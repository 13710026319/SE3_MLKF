% % =========================================================================
% % ESEKF、CUKF 与 CMLKF 对比测试脚本
% % 运行机制：加载指定的 4辆车和4个基站的 3D 数据进行仿真评估
% % 输出：控制台输出三算法联合 RMSE 评估表，作图显示各车欧氏定位误差对比
% % =========================================================================
clc; clear; close all;
%% 1. 加载路径与仿真数据
addpath(genpath('../Common'));
addpath(genpath('../Filter'));
addpath(genpath('../Data'));
% 导入指定的仿真数据集 (含4基站，提升X,Y方向可观测性)
data_file = 'E:\SE3_MLKF\Data\Trj_data_Veh4_Anc4_3D.mat';
if ~exist(data_file, 'file')
    % 若文件不存在，可备用查找相对路径
    data_file = '../Data/Trj_data_Veh4_Anc4_3D.mat';
    if ~exist(data_file, 'file')
        error('未检测到指定数据文件，请先运行 Data 下仿真生成脚本！');
    end
end
load(data_file); 
%% 2. 初始化滤波器 (同参数注入)
N_steps = length(trajectories.V1.Time_true);
dt_imu = 0.01;
% 标称状态初值生成
init_states = struct('X', {}, 'b', {});
for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    R0 = trajectories.(v_name).R_true(:, :, 1);
    v0 = [trajectories.(v_name).Vx_true(1); trajectories.(v_name).Vy_true(1); trajectories.(v_name).Vz_true(1)];
    p0 = [trajectories.(v_name).X_true(1); trajectories.(v_name).Y_true(1); trajectories.(v_name).Z_true(1)];
    
    X0 = eye(5);
    X0(1:3, 1:3) = R0;
    X0(1:3, 4)   = v0;
    X0(1:3, 5)   = p0;
    
    init_states(n).X = X0;
    init_states(n).b = zeros(6, 1);
end
% 协方差初值生成 [10]
P_n_init = diag([ ...
    (1*pi/180)^2 * ones(1, 3), ... % 姿态不确定
    (0.01)^2 * ones(1, 3), ...     % 速度不确定
    (0.01)^2 * ones(1, 3), ...     % 位置不确定
    (0.15)^2 * ones(1, 3), ...     % 加速度零偏初始范围 (扩充提升收敛速度)
    (0.015)^2 * ones(1, 3) ...     % 陀螺仪零偏初始范围
]);
init_P = kron(eye(Vehicle_num), P_n_init);
% 实例化三个滤波器，分别拷贝初值
filter_esekf = ESEKF(init_states, init_P, IMU_noise_params);
filter_cukf  = CUKF(init_states, init_P, IMU_noise_params);
filter_cmlkf = CMLKF(init_states, init_P, IMU_noise_params);
%% 3. 运行三滤波器混合仿真循环
pos_est_esekf = cell(Vehicle_num, 1);
pos_est_cukf  = cell(Vehicle_num, 1);
pos_est_cmlkf = cell(Vehicle_num, 1);
for n = 1:Vehicle_num
    pos_est_esekf{n} = zeros(N_steps, 3);
    pos_est_cukf{n}  = zeros(N_steps, 3);
    pos_est_cmlkf{n} = zeros(N_steps, 3);
    pos_est_esekf{n}(1, :) = init_states(n).X(1:3, 5)';
    pos_est_cukf{n}(1, :)  = init_states(n).X(1:3, 5)';
    pos_est_cmlkf{n}(1, :) = init_states(n).X(1:3, 5)';
end
fprintf('开始执行 ESEKF、CUKF 与 CMLKF 三算法对比运算...\n');
tic;
for k = 2:N_steps
    % 提取 IMU 观测
    imu_acc = zeros(Vehicle_num, 3);
    imu_gyro = zeros(Vehicle_num, 3);
    for n = 1:Vehicle_num
        v_name = sprintf('V%d', n);
        imu_acc(n, :) = trajectories.(v_name).IMU_acc_m(k, :);
        imu_gyro(n, :) = trajectories.(v_name).IMU_gyro_m(k, :);
    end
    
    % A. 预测更新步 (共用父类 manifold_predict 方法) [7]
    filter_esekf.propagate(imu_acc, imu_gyro, dt_imu);
    filter_cukf.propagate(imu_acc, imu_gyro, dt_imu);   % <--- 已补充 CUKF 预测
    filter_cmlkf.propagate(imu_acc, imu_gyro, dt_imu);
    
    % B. UWB 10Hz 观测纠正步
    if mod(k - 1, 10) == 0
        uwb_idx = (k - 1) / 10 + 1;
        
        anc_meas = zeros(Vehicle_num, Anchor_num);
        rel_meas = zeros(Vehicle_num, Vehicle_num);
        for n = 1:Vehicle_num
            v_name = sprintf('V%d', n);
            anc_meas(n, :) = trajectories.(v_name).UWB_Anchor(uwb_idx, 2:end);
            rel_meas(n, :) = trajectories.(v_name).UWB_Relative(uwb_idx, 2:end);
        end
        
        % 三种形式卡尔曼滤波更新
        filter_esekf.update(anchors, anc_meas, rel_meas, UWB_noise_params.sigma_anc, UWB_noise_params.sigma_rel);
        filter_cukf.update(anchors, anc_meas, rel_meas, UWB_noise_params.sigma_anc, UWB_noise_params.sigma_rel);    % <--- 已补充 CUKF 更新
        filter_cmlkf.update(anchors, anc_meas, rel_meas, UWB_noise_params.sigma_anc, UWB_noise_params.sigma_rel);
    end
    
    % C. 结果记录
    for n = 1:Vehicle_num
        pos_est_esekf{n}(k, :) = filter_esekf.states(n).X(1:3, 5)';
        pos_est_cukf{n}(k, :)  = filter_cukf.states(n).X(1:3, 5)';     % <--- 已补充 CUKF 状态记录
        pos_est_cmlkf{n}(k, :) = filter_cmlkf.states(n).X(1:3, 5)';
    end
end
toc;
fprintf('三滤波器运行完毕。\n');
%% 4. 位置误差与联合 RMSE 评估表格输出
pos_true = cell(Vehicle_num, 1);
for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    pos_true{n} = [trajectories.(v_name).X_true, trajectories.(v_name).Y_true, trajectories.(v_name).Z_true];
end
% 分别计算误差
[errors_esekf, rmse_esekf] = calculate_position_errors(pos_est_esekf, pos_true);
[errors_cukf, rmse_cukf]   = calculate_position_errors(pos_est_cukf, pos_true);     % <--- 已补充 CUKF 误差计算
[errors_cmlkf, rmse_cmlkf] = calculate_position_errors(pos_est_cmlkf, pos_true);
% 组装表格数据 [2]（行数扩充为 Vehicle_num * 3）
RowNames = cell(Vehicle_num * 3, 1);
X_RMSE = zeros(Vehicle_num * 3, 1);
Y_RMSE = zeros(Vehicle_num * 3, 1);
Z_RMSE = zeros(Vehicle_num * 3, 1);
Euc_RMSE = zeros(Vehicle_num * 3, 1);
for n = 1:Vehicle_num
    % ESEKF 行
    RowNames{3*n - 2} = sprintf('V%d_ESEKF', n);
    X_RMSE(3*n - 2)   = rmse_esekf(n).axis_rmse(1);
    Y_RMSE(3*n - 2)   = rmse_esekf(n).axis_rmse(2);
    Z_RMSE(3*n - 2)   = rmse_esekf(n).axis_rmse(3);
    Euc_RMSE(3*n - 2) = rmse_esekf(n).euc_rmse;
    
    % CUKF 行
    RowNames{3*n - 1} = sprintf('V%d_CUKF', n);
    X_RMSE(3*n - 1)   = rmse_cukf(n).axis_rmse(1);
    Y_RMSE(3*n - 1)   = rmse_cukf(n).axis_rmse(2);
    Z_RMSE(3*n - 1)   = rmse_cukf(n).axis_rmse(3);
    Euc_RMSE(3*n - 1) = rmse_cukf(n).euc_rmse;
    
    % CMLKF 行
    RowNames{3*n}     = sprintf('V%d_CMLKF', n);
    X_RMSE(3*n)       = rmse_cmlkf(n).axis_rmse(1);
    Y_RMSE(3*n)       = rmse_cmlkf(n).axis_rmse(2);
    Z_RMSE(3*n)       = rmse_cmlkf(n).axis_rmse(3);
    Euc_RMSE(3*n)     = rmse_cmlkf(n).euc_rmse;
end
% 打印表格
comp_table = table(X_RMSE, Y_RMSE, Z_RMSE, Euc_RMSE, 'RowNames', RowNames);
fprintf('\n======================= ESEKF VS CUKF VS CMLKF 联合评估表 =======================\n');
disp(comp_table);
fprintf('=================================================================================\n');
%% 5. 欧氏定位误差曲线对比图绘制 (无三轴误差)
time_arr = trajectories.V1.IMU_Time;
figure('Name', 'Euclidean Position Errors: ESEKF VS CUKF VS CMLKF', 'Position', [150, 150, 1000, 800]);
for n = 1:Vehicle_num
    subplot(2, 2, n);
    hold on; grid on;
    
    % 绘制对比线（新增绿色点划线代表 CUKF）
    plot(time_arr, errors_esekf(n).euc_err, 'b-', 'LineWidth', 1.2, 'DisplayName', 'ESEKF');
    plot(time_arr, errors_cukf(n).euc_err, 'g-.', 'LineWidth', 1.3, 'DisplayName', 'CUKF');
    plot(time_arr, errors_cmlkf(n).euc_err, 'r--', 'LineWidth', 1.5, 'DisplayName', 'CMLKF');
    
    title(sprintf('Vehicle %d Euclidean Error', n));
    xlabel('Time (s)');
    ylabel('Position Error (m)');
    xlim([0, time_arr(end)]);
    ylim([0, 1.5]); % 保持误差坐标视觉统一
    legend('Location', 'northeast');
    hold off;
end