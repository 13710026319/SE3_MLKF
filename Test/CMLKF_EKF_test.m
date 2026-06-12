% =========================================================================
% Filter_compare.m (极速 3D 无阻碍运行版)
% 核心对比：CMLKF_15D (100Hz高频MLKF更新) VS EKF_9D (10Hz经典更新)
% =========================================================================

clc; clear; close all;

%% 1. 加载路径与仿真数据
addpath(genpath('../Common'));
addpath(genpath('../Filter'));
addpath(genpath('../Data'));

data_file = 'E:\SE3_MLKF\Data\Trj_data_Veh4_Anc4_3D.mat';
if ~exist(data_file, 'file')
    data_file = '../Data/Trj_data_Veh4_Anc4_3D.mat'; 
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

%% 3. 初始化两大滤波器状态与协方差矩阵
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
    (0.01)^2 * ones(1, 3), ...     
    (0.01)^2 * ones(1, 3), ...     
    (0.05)^2 * ones(1, 3), ...     
    (1*pi/180)^2 * ones(1, 3), ... 
    (0.005)^2 * ones(1, 3) ...     
]);
init_P_15d = kron(eye(Vehicle_num), P_n_init_15d);

Q_sigmas_15d.sig_wp      = 0;  
Q_sigmas_15d.sig_wv      = 0;  
Q_sigmas_15d.sig_wa      = 0.0005; 
Q_sigmas_15d.sig_wR      = 0.00005;  
Q_sigmas_15d.sig_womega = 0.00005;


% B. 为 9 维经典算法组 (EKF) 初始化标称状态与协方差 [1.1.5]
init_states_9d = struct('p', {}, 'v', {}, 'R', {});
for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    veh = trajectories.(v_name);
    init_states_9d(n).p = [veh.X_true(1); veh.Y_true(1); veh.Z_true(1)];
    init_states_9d(n).v = [veh.Vx_true(1); veh.Vy_true(1); veh.Vz_true(1)];
    init_states_9d(n).R = veh.R_true(:, :, 1);
end

P_n_init_9d = diag([ ...
    (0.01)^2 * ones(1, 3), ...     
    (0.01)^2 * ones(1, 3), ...     
    (1*pi/180)^2 * ones(1, 3) ...  
]);
init_P_9d = kron(eye(Vehicle_num), P_n_init_9d);

% C. 实例化两个滤波器类
filter_cmlkf = CMLKF(init_states_15d, init_P_15d, Q_sigmas_15d);
filter_ekf   = EKF(init_states_9d, init_P_9d); 

%% 4. 执行多滤波器并行仿真循环
pos_est_cmlkf = cell(Vehicle_num, 1);
pos_est_ekf   = cell(Vehicle_num, 1);
for n = 1:Vehicle_num
    pos_est_cmlkf{n} = zeros(N_steps, 3);
    pos_est_ekf{n}   = zeros(N_steps, 3);
    pos_est_cmlkf{n}(1, :) = init_states_15d(n).p';
    pos_est_ekf{n}(1, :)   = init_states_9d(n).p';
end

fprintf('开始执行 CMLKF 与 EKF 协同估计运行...\n');
tic;
for k = 2:N_steps
    % 提取当前 IMU 真实零偏并预校正 [1]
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
    
    % A. 预测步：状态时间传播
    filter_cmlkf.propagate(dt_imu);
    filter_ekf.propagate(imu_acc, imu_gyro, dt_imu); % EKF 9D 以 IMU 作为控制输入 [1.2.9]
    
    % B. 观测更新步 [8]
    if mod(k - 1, 10) == 0
        % --- UWB测距联合更新周期 (10Hz) ---
        uwb_idx = (k - 1) / 10 + 1;
        
        anc_meas = zeros(Vehicle_num, Anchor_num);
        rel_meas = zeros(Vehicle_num, Vehicle_num);
        for n = 1:Vehicle_num
            v_name = sprintf('V%d', n);
            veh = trajectories.(v_name);
            anc_meas(n, :) = veh.UWB_Anchor(uwb_idx, 2:end);
            rel_meas(n, :) = veh.UWB_Relative(uwb_idx, 2:end);
        end
        
        % CMLKF：15D 联合最大似然估计更新步 (参数已正确对齐 10 个接口) [4]
        filter_cmlkf.update(imu_acc, imu_gyro, anchors, anc_meas, rel_meas, ...
                            IMU_noise_params.sigma_na, IMU_noise_params.sigma_nw, ...
                            UWB_noise_params.sigma_anc, UWB_noise_params.sigma_rel);
                      
        % EKF (9D)：IMU作为输入，观测步仅依靠 UWB 更新 [1.1.5, 1.2.9]
        filter_ekf.update(anchors, anc_meas, rel_meas, ...
                          UWB_noise_params.sigma_anc, UWB_noise_params.sigma_rel);
    else
        % --- 高频 IMU-only 最大似然更新步 (其余 90Hz 步) ---
        filter_cmlkf.update_imu_only(imu_acc, imu_gyro, ...
                                     IMU_noise_params.sigma_na, IMU_noise_params.sigma_nw);
    end
    
    % C. 结果记录
    for n = 1:Vehicle_num
        pos_est_cmlkf{n}(k, :) = filter_cmlkf.states(n).p';
        pos_est_ekf{n}(k, :)   = filter_ekf.states(n).p';
    end
end
toc;
fprintf('双算法仿真对比计算完毕。\n');

%% 5. 位置误差计算与控制台表格联合输出
pos_true = cell(Vehicle_num, 1);
for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    pos_true{n} = [trajectories.(v_name).X_true, ...
                   trajectories.(v_name).Y_true, ...
                   trajectories.(v_name).Z_true];
end

[errors_cmlkf, rmse_cmlkf] = calculate_position_errors(pos_est_cmlkf, pos_true);
[errors_ekf, rmse_ekf]     = calculate_position_errors(pos_est_ekf, pos_true);

% 组装比较表 [2]
RowNames = cell(Vehicle_num * 2, 1);
X_RMSE = zeros(Vehicle_num * 2, 1);
Y_RMSE = zeros(Vehicle_num * 2, 1);
Z_RMSE = zeros(Vehicle_num * 2, 1);
Euc_RMSE = zeros(Vehicle_num * 2, 1);

for n = 1:Vehicle_num
    % 1. EKF (9D)
    idx_ekf = 2*n - 1;
    RowNames{idx_ekf} = sprintf('V%d_EKF_9D(IMU_Input)', n);
    X_RMSE(idx_ekf)   = rmse_ekf(n).axis_rmse(1);
    Y_RMSE(idx_ekf)   = rmse_ekf(n).axis_rmse(2);
    Z_RMSE(idx_ekf)   = rmse_ekf(n).axis_rmse(3);
    Euc_RMSE(idx_ekf) = rmse_ekf(n).euc_rmse;
    
    % 2. CMLKF (15D，重构优化版)
    idx_new_ml = 2*n;
    RowNames{idx_new_ml}   = sprintf('V%d_CMLKF_15D(IMU_Obs_100Hz_MLKF)', n);
    X_RMSE(idx_new_ml)     = rmse_cmlkf(n).axis_rmse(1);
    Y_RMSE(idx_new_ml)     = rmse_cmlkf(n).axis_rmse(2);
    Z_RMSE(idx_new_ml)     = rmse_cmlkf(n).axis_rmse(3);
    Euc_RMSE(idx_new_ml)   = rmse_cmlkf(n).euc_rmse;
end

rmse_comparison_table = table(X_RMSE, Y_RMSE, Z_RMSE, Euc_RMSE, 'RowNames', RowNames);
fprintf('\n========================== CMLKF (流形高频版) VS EKF (离散噪声版) 性能评估对比表 ==========================\n');
disp(rmse_comparison_table);
fprintf('===================================================================================================\n');

%% 6. 绘图对比
time_arr = trajectories.V1.IMU_Time;
figure('Name', 'Euclidean Position Errors: CMLKF VS EKF', 'Position', [100, 100, 1100, 800]);

for n = 1:Vehicle_num
    subplot(2, 2, n);
    hold on; grid on;
    
    % 1. 经典 9维 EKF (蓝色实线) [1.1.5]
    plot(time_arr, errors_ekf(n).euc_err, 'b-', 'LineWidth', 1.2, 'DisplayName', 'Centralized EKF (9D, IMU as Input)');
    % 2. 优化版 CMLKF (红色虚线 - 现已应用100Hz最大似然流形优化更新 + J_t=I) [4, 7]
    plot(time_arr, errors_cmlkf(n).euc_err, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Centralized MLKF (15D, MLKF-IMU, 100Hz, J_t=I)');
    
    title(sprintf('Vehicle %d Euclidean Error Comparison', n));
    xlabel('Time (s)'); ylabel('Error (m)');
    xlim([0, time_arr(end)]);
    ylim([0, 1.5]); 
    legend('Location', 'northeast');
    hold off;
end