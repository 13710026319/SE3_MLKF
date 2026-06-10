% =========================================================================
% CMLKF_VS_EKF.m
% 新版3D多车(无人机)流形算法对比测试脚本：CMLKF VS 经典集中式 EKF (ESEKF)
% 运行机制：双滤波器在相同初值、噪声和IMU/UWB同步数据下并排运算 [3, 4, 8]
% 输出展示：控制台输出联合 RMSE 评估表，图形窗口输出各车欧氏定位误差对比曲线 [2]
% =========================================================================

clc; clear; close all;

%% 1. 加载路径与仿真数据
addpath(genpath('../Common'));
addpath(genpath('../Filter'));
addpath(genpath('../Data'));

% 导入平滑转向后的3D基站不共面数据集
data_file = 'E:\SE3_MLKF\Data\Trj_data_Veh4_Anc4_3D.mat';
if ~exist(data_file, 'file')
    data_file = '../Data/Trj_data_Veh4_Anc4_3D.mat'; % 相对路径备用
    if ~exist(data_file, 'file')
        error('未检测到数据集文件，请确保已运行 generate_data_3D.m 生成数据！');
    end
end
load(data_file); % 加载 trajectories, anchors, IMU_noise_params, UWB_noise_params, Vehicle_num, Anchor_num

dt_imu = 0.01; % 100Hz 仿真步长 \tau [1]

%% 2. 状态真值重建与偏置已知设定
% 离线计算标称状态中的加速度与角速度真值，保证初始化一致性 [1, 3]
for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    veh = trajectories.(v_name);
    N_steps = length(veh.Time_true);
    
    % A. 重建 3D 真实加速度 a_t^i
    v_true_matrix = [veh.Vx_true, veh.Vy_true, veh.Vz_true];
    a_true_matrix = zeros(N_steps, 3);
    a_true_matrix(:, 1) = gradient(v_true_matrix(:, 1), dt_imu);
    a_true_matrix(:, 2) = gradient(v_true_matrix(:, 2), dt_imu);
    a_true_matrix(:, 3) = gradient(v_true_matrix(:, 3), dt_imu);
    trajectories.(v_name).a_true = a_true_matrix;
    
    % B. 重建 3D 真实角速度 \omega_t^i
    theta_unwrapped = unwrap(veh.Theta_true);
    wz_true = gradient(theta_unwrapped, dt_imu);
    omega_true_matrix = [zeros(N_steps, 2), wz_true];
    trajectories.(v_name).omega_true = omega_true_matrix;
end

%% 3. 初始化 CMLKF 与 EKF 滤波器 (使用完全相同的参数)
init_states = struct('p', {}, 'v', {}, 'a', {}, 'R', {}, 'omega', {});
for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    veh = trajectories.(v_name);
    
    % 提取 3D 动力学状态初值 [1]
    init_states(n).p     = [veh.X_true(1); veh.Y_true(1); veh.Z_true(1)];
    init_states(n).v     = [veh.Vx_true(1); veh.Vy_true(1); veh.Vz_true(1)];
    init_states(n).a     = veh.a_true(1, :)';
    init_states(n).R     = veh.R_true(:, :, 1);
    init_states(n).omega = veh.omega_true(1, :)';
end

% 构造单车误差状态初始协方差 P_n (15x15) [2, 3]
P_n_init = diag([ ...
    (0.01)^2 * ones(1, 3), ...     % 位置不确定度
    (0.01)^2 * ones(1, 3), ...     % 速度不确定度
    (0.05)^2 * ones(1, 3), ...     % 加速度不确定度
    (1*pi/180)^2 * ones(1, 3), ... % 姿态不确定度 (1 deg)
    (0.005)^2 * ones(1, 3) ...     % 角速度不确定度
]);
init_P = kron(eye(Vehicle_num), P_n_init);

% 过程噪声标准差设置 [1]
Q_sigmas.sig_wp     = 0.001;  
Q_sigmas.sig_wv     = 0.005;  
Q_sigmas.sig_wa     = 0.05;   
Q_sigmas.sig_wR     = 0.001;  
Q_sigmas.sig_womega = 0.005;  

% 独立实例化两个滤波器，拷贝完全一致的初值与噪声配置
filter_cmlkf = CMLKF(init_states, init_P, Q_sigmas);
filter_ekf   = EKF(init_states, init_P, Q_sigmas);

%% 4. 执行双滤波器并行仿真循环
pos_est_cmlkf = cell(Vehicle_num, 1);
pos_est_ekf   = cell(Vehicle_num, 1);
for n = 1:Vehicle_num
    pos_est_cmlkf{n} = zeros(N_steps, 3);
    pos_est_ekf{n}   = zeros(N_steps, 3);
    pos_est_cmlkf{n}(1, :) = init_states(n).p';
    pos_est_ekf{n}(1, :)   = init_states(n).p';
end

fprintf('开始执行双滤波器 (CMLKF VS EKF) 状态估计对比...\n');
tic;
for k = 2:N_steps
    % A. 预测步：纯离散动力学状态时间传播与协方差前向传播 [3, 4]
    filter_cmlkf.propagate(dt_imu);
    filter_ekf.propagate(dt_imu);
    
    % B. UWB 测距与 IMU 联合更新步 (10Hz 严格同步周期) [8]
    if mod(k - 1, 10) == 0
        uwb_idx = (k - 1) / 10 + 1;
        
        imu_acc = zeros(Vehicle_num, 3);
        imu_gyro = zeros(Vehicle_num, 3);
        anc_meas = zeros(Vehicle_num, Anchor_num);
        rel_meas = zeros(Vehicle_num, Vehicle_num);
        
        for n = 1:Vehicle_num
            v_name = sprintf('V%d', n);
            veh = trajectories.(v_name);
            
            % 读取真实的 IMU 零偏进行已知常数扣除 [1]
            ba_true = veh.IMU_bias_a_true(k, :)';
            bw_true = veh.IMU_bias_w_true(k, :)';
            
            % 预修正 IMU 观测： \tilde{a} 和 \tilde{\omega} (文档 Eq 7, 8) [1]
            imu_acc(n, :)  = (veh.IMU_acc_m(k, :)'  - ba_true)';
            imu_gyro(n, :) = (veh.IMU_gyro_m(k, :)' - bw_true)';
            
            % 提取当前周期 UWB 观测值 [8]
            anc_meas(n, :) = veh.UWB_Anchor(uwb_idx, 2:end);
            rel_meas(n, :) = veh.UWB_Relative(uwb_idx, 2:end);
        end
        
        % CMLKF 观测更新 (通过测量流形非线性 Gauss-Newton 迭代求解) [4]
        filter_cmlkf.update(imu_acc, imu_gyro, anchors, anc_meas, rel_meas, ...
                            IMU_noise_params.sigma_na, IMU_noise_params.sigma_nw, ...
                            UWB_noise_params.sigma_anc, UWB_noise_params.sigma_rel);
                        
        % EKF 观测更新 (直接在先验工作点一阶线性化更新) [4]
        filter_ekf.update(imu_acc, imu_gyro, anchors, anc_meas, rel_meas, ...
                          IMU_noise_params.sigma_na, IMU_noise_params.sigma_nw, ...
                          UWB_noise_params.sigma_anc, UWB_noise_params.sigma_rel);
    end
    
    % C. 结果提取与记录
    for n = 1:Vehicle_num
        pos_est_cmlkf{n}(k, :) = filter_cmlkf.states(n).p';
        pos_est_ekf{n}(k, :)   = filter_ekf.states(n).p';
    end
end
toc;
fprintf('双算法仿真对比计算完毕。\n');

%% 5. 误差计算与 RMSE 表格联合输出
pos_true = cell(Vehicle_num, 1);
for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    pos_true{n} = [trajectories.(v_name).X_true, ...
                   trajectories.(v_name).Y_true, ...
                   trajectories.(v_name).Z_true];
end

% 分别分析两个算法的位置误差
[errors_cmlkf, rmse_cmlkf] = calculate_position_errors(pos_est_cmlkf, pos_true);
[errors_ekf, rmse_ekf]     = calculate_position_errors(pos_est_ekf, pos_true);

% 组装对比表格结构 [2]
RowNames = cell(Vehicle_num * 2, 1);
X_RMSE = zeros(Vehicle_num * 2, 1);
Y_RMSE = zeros(Vehicle_num * 2, 1);
Z_RMSE = zeros(Vehicle_num * 2, 1);
Euc_RMSE = zeros(Vehicle_num * 2, 1);

for n = 1:Vehicle_num
    % EKF 算法项
    RowNames{2*n - 1} = sprintf('V%d_EKF', n);
    X_RMSE(2*n - 1)   = rmse_ekf(n).axis_rmse(1);
    Y_RMSE(2*n - 1)   = rmse_ekf(n).axis_rmse(2);
    Z_RMSE(2*n - 1)   = rmse_ekf(n).axis_rmse(3);
    Euc_RMSE(2*n - 1) = rmse_ekf(n).euc_rmse;
    
    % CMLKF 算法项
    RowNames{2*n}     = sprintf('V%d_CMLKF', n);
    X_RMSE(2*n)       = rmse_cmlkf(n).axis_rmse(1);
    Y_RMSE(2*n)       = rmse_cmlkf(n).axis_rmse(2);
    Z_RMSE(2*n)       = rmse_cmlkf(n).axis_rmse(3);
    Euc_RMSE(2*n)     = rmse_cmlkf(n).euc_rmse;
end

rmse_comp_table = table(X_RMSE, Y_RMSE, Z_RMSE, Euc_RMSE, 'RowNames', RowNames);
fprintf('\n======================== CMLKF VS EKF 定位性能评估表 ========================\n');
disp(rmse_comp_table);
fprintf('=============================================================================\n');

%% 6. 绘图 (仅对比欧氏定位误差曲线)
time_arr = trajectories.V1.IMU_Time;
figure('Name', 'Euclidean Position Errors: CMLKF VS EKF', 'Position', [150, 150, 1000, 800]);

for n = 1:Vehicle_num
    subplot(2, 2, n);
    hold on; grid on;
    
    % 绘制 EKF 欧氏定位误差线 (蓝色实线)
    plot(time_arr, errors_ekf(n).euc_err, 'b-', 'LineWidth', 1.2, 'DisplayName', 'Centralized EKF');
    % 绘制 CMLKF 欧氏定位误差线 (红色虚线)
    plot(time_arr, errors_cmlkf(n).euc_err, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Centralized MLKF');
    
    title(sprintf('Vehicle %d Euclidean Position Error', n));
    xlabel('Time (s)'); ylabel('Error (m)');
    xlim([0, time_arr(end)]);
    ylim([0, 1.5]); % 坐标轴统一限制，方便肉眼直观分析
    legend('Location', 'northeast');
    hold off;
end