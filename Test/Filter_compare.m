% =========================================================================
% Filter_compare.m
% 双算法并排对比测试脚本：CMLKF VS 经典 ESKF (IMU作为输入)
% 运行机制：加载指定的 3D 数据集，在完全相同的 IMU 偏置和 UWB 测距下并排运算
% 改进重点：为 CMLKF 开启高频 IMU-only 卡尔曼更新，充分挖掘 100Hz 的 IMU 数据性能 [1, 3]
% =========================================================================

clc; clear; close all;

%% 1. 加载路径与仿真数据
addpath(genpath('../Common'));
addpath(genpath('../Filter'));
addpath(genpath('../Data'));

% 导入指定的仿真数据集 (4基站，高度不共面) [1.1.9]
data_file = 'E:\SE3_MLKF\Data\Trj_data_Veh4_Anc4_3D.mat';
if ~exist(data_file, 'file')
    data_file = '../Data/Trj_data_Veh4_Anc4_3D.mat'; % 相对路径备用
    if ~exist(data_file, 'file')
        error('未检测到指定数据文件，请先运行 Data 下数据生成脚本！');
    end
end
load(data_file); % 加载获得 trajectories, anchors, IMU_noise_params, UWB_noise_params, Vehicle_num, Anchor_num

dt_imu = 0.01; % 100Hz 仿真步长 \tau [1]

%% 2. 状态真值重建与偏置已知设定
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

%% 3. 初始化两大滤波器状态与协方差矩阵 (严格对齐物理假设)

% A. 为 15 维算法组 (CMLKF) 初始化标称状态与协方差 [1, 2, 3]
init_states_15d = struct('p', {}, 'v', {}, 'a', {}, 'R', {}, 'omega', {});
for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    veh = trajectories.(v_name);
    init_states_15d(n).p     = [veh.X_true(1); veh.Y_true(1); veh.Z_true(1)];
    init_states_15d(n).v     = [veh.Vx_true(1); veh.Vy_true(1); veh.Vz_true(1)];
    init_states_15d(n).a     = veh.a_true(1, :)';
    init_states_15d(n).R     = veh.R_true(:, :, 1);
    init_states_15d(n).omega = veh.omega_true(1, :)';
end

P_n_init_15d = diag([ ...
(0.02)^2 * ones(1, 3), ... % 位置不确定度
(0.02)^2 * ones(1, 3), ... % 速度不确定度
(0.2)^2 * ones(1, 3), ... % 加速度不确定度
(2*pi/180)^2 * ones(1, 3), ... % 姿态不确定度
(0.02)^2 * ones(1, 3) ... % 角速度不确定度
]);

init_P_15d = kron(eye(Vehicle_num), P_n_init_15d);

% --- 核心修改2：降低 az 的过程噪声 ---
Q_sigmas_15d.sig_wp      = 0.005;  
Q_sigmas_15d.sig_wv      = 0.015;  
Q_sigmas_15d.sig_wa      = 0.025; 
Q_sigmas_15d.sig_wR      = 0.005;  
Q_sigmas_15d.sig_womega = 0.0025;

% B. 为 9 维经典算法组 (ESKF) 初始化标称状态与协方差 [1.1.5]
init_states_9d = struct('p', {}, 'v', {}, 'R', {});
for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    veh = trajectories.(v_name);
    init_states_9d(n).p = [veh.X_true(1); veh.Y_true(1); veh.Z_true(1)];
    init_states_9d(n).v = [veh.Vx_true(1); veh.Vy_true(1); veh.Vz_true(1)];
    init_states_9d(n).R = veh.R_true(:, :, 1);
end

P_n_init_9d = diag([ ...
    (0.01)^2 * ones(1, 3), ...     % 位置不确定度
    (0.01)^2 * ones(1, 3), ...     % 速度不确定度
    (1*pi/180)^2 * ones(1, 3) ...  % 姿态不确定度
]);
init_P_9d = kron(eye(Vehicle_num), P_n_init_9d);

% C. 实例化两个滤波器类
filter_cmlkf = CMLKF(init_states_15d, init_P_15d, Q_sigmas_15d);
filter_eskf  = ESKF(init_states_9d, init_P_9d, IMU_noise_params);


%% 4. 执行多滤波器并行仿真循环
pos_est_cmlkf = cell(Vehicle_num, 1);
pos_est_eskf  = cell(Vehicle_num, 1);
for n = 1:Vehicle_num
    pos_est_cmlkf{n} = zeros(N_steps, 3);
    pos_est_eskf{n}  = zeros(N_steps, 3);
    pos_est_cmlkf{n}(1, :) = init_states_15d(n).p';
    pos_est_eskf{n}(1, :)  = init_states_9d(n).p';
end

fprintf('开始执行双算法并排对比仿真运算...\n');
tic;
for k = 2:N_steps
    % 提取当前 IMU 真实零偏，对 IMU 读数进行高精度离线预修正 [1]
    imu_acc = zeros(Vehicle_num, 3);
    imu_gyro = zeros(Vehicle_num, 3);
    for n = 1:Vehicle_num
        v_name = sprintf('V%d', n);
        veh = trajectories.(v_name);
        ba_true = veh.IMU_bias_a_true(k, :)';
        bw_true = veh.IMU_bias_w_true(k, :)';
        imu_acc(n, :)  = (veh.IMU_acc_m(k, :)'  - ba_true)';
        imu_gyro(n, :) = (veh.IMU_gyro_m(k, :)' - bw_true)';
    end
    
    % A. 预测步：状态传播 [3]
    % 15 维 CMLKF 依靠系统动力学进行时间传播 [3]
    filter_cmlkf.propagate(dt_imu);
    % 9 维 ESKF 严格以预修正 IMU 数据为控制输入驱动状态时间传播 [1.2.9]
    filter_eskf.propagate(imu_acc, imu_gyro, dt_imu);
    
    % B. 观测更新步 [8]
    if mod(k - 1, 10) == 0
        % --- UWB测距更新周期 (10Hz) ---
        uwb_idx = (k - 1) / 10 + 1;
        
        anc_meas = zeros(Vehicle_num, Anchor_num);
        rel_meas = zeros(Vehicle_num, Vehicle_num);
        for n = 1:Vehicle_num
            v_name = sprintf('V%d', n);
            veh = trajectories.(v_name);
            anc_meas(n, :) = veh.UWB_Anchor(uwb_idx, 2:end);
            rel_meas(n, :) = veh.UWB_Relative(uwb_idx, 2:end);
        end
        
        % CMLKF：测量流形联合最大似然估计非线性投影更新 [4]
        filter_cmlkf.update(imu_acc, imu_gyro, anchors, anc_meas, rel_meas, ...
                            IMU_noise_params.sigma_na, IMU_noise_params.sigma_nw, ...
                            UWB_noise_params.sigma_anc, UWB_noise_params.sigma_rel);
                        
        % ESKF (9D)：IMU已做系统输入，观测步仅依靠 UWB 测距信息更新 [1.1.5, 1.2.9]
        filter_eskf.update(anchors, anc_meas, rel_meas, ...
                           UWB_noise_params.sigma_anc, UWB_noise_params.sigma_rel);
    else
        
        filter_cmlkf.update_imu_only(imu_acc, imu_gyro, ...
                                     IMU_noise_params.sigma_na, IMU_noise_params.sigma_nw);
    end
    
    % C. 结果记录
    for n = 1:Vehicle_num
        pos_est_cmlkf{n}(k, :) = filter_cmlkf.states(n).p';
        pos_est_eskf{n}(k, :)  = filter_eskf.states(n).p';
    end
end
toc;
fprintf('双算法并排估计计算完毕。\n');

%% 5. 位置误差计算与控制台表格联合输出
pos_true = cell(Vehicle_num, 1);
for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    pos_true{n} = [trajectories.(v_name).X_true, ...
                   trajectories.(v_name).Y_true, ...
                   trajectories.(v_name).Z_true];
end

[errors_cmlkf, rmse_cmlkf] = calculate_position_errors(pos_est_cmlkf, pos_true);
[errors_eskf, rmse_eskf]   = calculate_position_errors(pos_est_eskf, pos_true);

% 组装双算法联合 RMSE 比较表 [2]
RowNames = cell(Vehicle_num * 2, 1);
X_RMSE = zeros(Vehicle_num * 2, 1);
Y_RMSE = zeros(Vehicle_num * 2, 1);
Z_RMSE = zeros(Vehicle_num * 2, 1);
Euc_RMSE = zeros(Vehicle_num * 2, 1);

for n = 1:Vehicle_num
    % 1. ESKF (9D, 惯导输入架构)
    idx_eskf = 2*n - 1;
    RowNames{idx_eskf} = sprintf('V%d_ESKF_9D(IMU_Input)', n);
    X_RMSE(idx_eskf)   = rmse_eskf(n).axis_rmse(1);
    Y_RMSE(idx_eskf)   = rmse_eskf(n).axis_rmse(2);
    Z_RMSE(idx_eskf)   = rmse_eskf(n).axis_rmse(3);
    Euc_RMSE(idx_eskf) = rmse_eskf(n).euc_rmse;
    
    % 2. CMLKF (15D, 开启100Hz高频IMU-only更新)
    idx_ml = 2*n;
    RowNames{idx_ml}   = sprintf('V%d_CMLKF_15D(IMU_Obs_100Hz)', n);
    X_RMSE(idx_ml)     = rmse_cmlkf(n).axis_rmse(1);
    Y_RMSE(idx_ml)     = rmse_cmlkf(n).axis_rmse(2);
    Z_RMSE(idx_ml)     = rmse_cmlkf(n).axis_rmse(3);
    Euc_RMSE(idx_ml)   = rmse_cmlkf(n).euc_rmse;
end

rmse_comparison_table = table(X_RMSE, Y_RMSE, Z_RMSE, Euc_RMSE, 'RowNames', RowNames);
fprintf('\n========================== CMLKF VS ESKF 性能联合评估表 ==========================\n');
disp(rmse_comparison_table);
fprintf('=========================================================================================\n');

%% 6. 绘图对比 (3D 欧氏定位误差对比曲线)
time_arr = trajectories.V1.IMU_Time;
figure('Name', 'Euclidean Position Errors: ESKF VS CMLKF', 'Position', [100, 100, 1100, 800]);

for n = 1:Vehicle_num
    subplot(2, 2, n);
    hold on; grid on;
    
    % 1. 经典 9维 ESKF (绿色点划线) [1.1.5]
    plot(time_arr, errors_eskf(n).euc_err, 'g-.', 'LineWidth', 1.3, 'DisplayName', 'Centralized ESKF (9D, IMU as Input)');
    % 2. 15维 CMLKF (红色虚线 - 现已应用100Hz高频更新)
    plot(time_arr, errors_cmlkf(n).euc_err, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Centralized MLKF (15D, IMU as Obs, 100Hz)');
    
    title(sprintf('Vehicle %d Euclidean Error Comparison', n));
    xlabel('Time (s)'); ylabel('Error (m)');
    xlim([0, time_arr(end)]);
    ylim([0, 1.5]); 
    legend('Location', 'northeast');
    hold off;
end