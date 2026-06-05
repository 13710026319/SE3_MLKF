% =========================================================================
% ESEKF 算法运行测试脚本
% 运行机制：IMU 100Hz 传播预测，UWB 10Hz 触发观测更新步
% 输出展示：控制台输出 RMSE，图形输出 2x2 各车误差曲线图
% =========================================================================

clc; clear; close all;

%% 1. 加载路径与仿真数据
% 使用相对路径将 Common、Filter 和 Data 添加进 MATLAB 工作环境
addpath(genpath('../Common'));
addpath(genpath('../Filter'));
addpath(genpath('../Data'));

% 数据文件路径 (根据生成文件名加载)
data_file = 'E:\SE3_MLKF\Data\Trj_data_Veh4_Anc4_3D.mat';
if ~exist(data_file, 'file')
    error('未检测到数据文件，请先运行 Data 文件夹下的仿真生成脚本！');
end
load(data_file); % 加载获得 trajectories, anchors, IMU_noise_params, UWB_noise_params, Vehicle_num, Anchor_num

%% 2. 初始化 ESEKF 滤波器
% 获取总时间步数
N_steps = length(trajectories.V1.Time_true);
dt_imu = 0.01; % 100Hz

% 构造滤波器初始状态结构体数组
init_states = struct('X', {}, 'b', {});
for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    
    % 从首帧数据中获取真值作为滤波器起点
    R0 = trajectories.(v_name).R_true(:, :, 1);
    v0 = [trajectories.(v_name).Vx_true(1); 
          trajectories.(v_name).Vy_true(1); 
          trajectories.(v_name).Vz_true(1)];
    p0 = [trajectories.(v_name).X_true(1); 
          trajectories.(v_name).Y_true(1); 
          trajectories.(v_name).Z_true(1)];
    
    % 构造 5x5 的 SE2(3) 矩阵
    X0 = eye(5);
    X0(1:3, 1:3) = R0;
    X0(1:3, 4)   = v0;
    X0(1:3, 5)   = p0;
    
    init_states(n).X = X0;
    init_states(n).b = zeros(6, 1); % 零偏初始假设为 0
end

% 构造单车误差状态初始协方差 P_n (15x15) [10]
P_n_init = diag([ ...
    (1*pi/180)^2 * ones(1, 3), ... % 姿态方差 (1 deg)
    (0.01)^2 * ones(1, 3), ...     % 速度方差 (0.01 m/s)
    (0.01)^2 * ones(1, 3), ...     % 位置方差 (0.01 m)
    (0.015)^2 * ones(1, 3), ...     % 加速度偏置方差
    (0.0015)^2 * ones(1, 3) ...     % 陀螺仪偏置方差
]);

% 拼装多车联合协方差矩阵 P (15N x 15N)
init_P = kron(eye(Vehicle_num), P_n_init);

% 实例化 ESEKF 滤波器
filter = ESEKF(init_states, init_P, IMU_noise_params);

%% 3. 滤波器主循环运行步
% 预分配估计位置存储器 (Cell 数组)
pos_est = cell(Vehicle_num, 1);
for n = 1:Vehicle_num
    pos_est{n} = zeros(N_steps, 3);
    pos_est{n}(1, :) = init_states(n).X(1:3, 5)';
end

fprintf('开始执行 ESEKF 滤波估计...\n');
tic;
for k = 2:N_steps
    % A. 提取所有车辆当前时间步的 3D IMU 数据
    imu_acc = zeros(Vehicle_num, 3);
    imu_gyro = zeros(Vehicle_num, 3);
    for n = 1:Vehicle_num
        v_name = sprintf('V%d', n);
        imu_acc(n, :) = trajectories.(v_name).IMU_acc_m(k, :);
        imu_gyro(n, :) = trajectories.(v_name).IMU_gyro_m(k, :);
    end
    
    % B. 状态和协方差传播 (IMU 100Hz 预测)
    filter.propagate(imu_acc, imu_gyro, dt_imu);
    
    % C. 判断是否为 10Hz 的 UWB 测量更新时刻
    % (由于 100Hz:10Hz 对齐，每 10 步触发一次 UWB 更新，对齐 1, 11, 21...)
    if mod(k - 1, 10) == 0
        uwb_idx = (k - 1) / 10 + 1;
        
        % 提取 UWB 基站观测与车间相对观测
        anc_meas = zeros(Vehicle_num, Anchor_num);
        rel_meas = zeros(Vehicle_num, Vehicle_num);
        for n = 1:Vehicle_num
            v_name = sprintf('V%d', n);
            anc_meas(n, :) = trajectories.(v_name).UWB_Anchor(uwb_idx, 2:end);
            rel_meas(n, :) = trajectories.(v_name).UWB_Relative(uwb_idx, 2:end);
        end
        
        % 执行 3D UWB 卡尔曼更新 [8]
        filter.update(anchors, anc_meas, rel_meas, ...
            UWB_noise_params.sigma_anc, UWB_noise_params.sigma_rel);
    end
    
    % D. 收集并记录本时间步的估计结果
    for n = 1:Vehicle_num
        pos_est{n}(k, :) = filter.states(n).X(1:3, 5)';
    end
end
toc;
fprintf('滤波计算完毕！\n');

%% 4. 位置误差与 RMSE 计算
% 提取真实位置用于评估对比
pos_true = cell(Vehicle_num, 1);
for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    pos_true{n} = [trajectories.(v_name).X_true, ...
                   trajectories.(v_name).Y_true, ...
                   trajectories.(v_name).Z_true];
end

% 调用误差分析工具
[errors, rmse] = calculate_position_errors(pos_est, pos_true);

% --- 修改为表格形式输出 RMSE ---
RowNames = cell(Vehicle_num, 1);
X_RMSE = zeros(Vehicle_num, 1);
Y_RMSE = zeros(Vehicle_num, 1);
Z_RMSE = zeros(Vehicle_num, 1);
Euc_RMSE = zeros(Vehicle_num, 1);

for n = 1:Vehicle_num
    RowNames{n} = sprintf('Vehicle_%d', n);
    X_RMSE(n)   = rmse(n).axis_rmse(1);
    Y_RMSE(n)   = rmse(n).axis_rmse(2);
    Z_RMSE(n)   = rmse(n).axis_rmse(3);
    Euc_RMSE(n) = rmse(n).euc_rmse;
end

% 组装 MATLAB 官方表格对象并打印 [2]
rmse_table = table(X_RMSE, Y_RMSE, Z_RMSE, Euc_RMSE, 'RowNames', RowNames);

fprintf('\n======================== ESEKF 定位 RMSE 评估表 ========================\n');
disp(rmse_table);
fprintf('========================================================================\n');

%% 5. 误差曲线绘图
time_arr = trajectories.V1.IMU_Time;

figure('Name', 'ESEKF Position Errors', 'Position', [150, 150, 1000, 800]);
for n = 1:Vehicle_num
    subplot(2, 2, n);
    hold on; grid on;
    
    % 绘制 X, Y, Z 三轴及欧氏距离误差曲线
    plot(time_arr, errors(n).axis_err(:, 1), 'r-', 'LineWidth', 1.2, 'DisplayName', 'X Error');
    plot(time_arr, errors(n).axis_err(:, 2), 'g-', 'LineWidth', 1.2, 'DisplayName', 'Y Error');
    plot(time_arr, errors(n).axis_err(:, 3), 'b-', 'LineWidth', 1.2, 'DisplayName', 'Z Error');
    plot(time_arr, errors(n).euc_err, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Euclidean Error');
    
    title(sprintf('Vehicle %d Position Error', n));
    xlabel('Time (s)');
    ylabel('Position Error (m)');
    xlim([0, time_arr(end)]);
    ylim([-1.5, 1.5]); % 根据噪声情况和收敛后范围，此范围通常利于观察收敛细节
    legend('Location', 'northeast');
    hold off;
end