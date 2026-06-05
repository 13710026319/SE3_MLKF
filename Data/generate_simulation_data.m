% =========================================================================
% SE2(3)-MLKF 算法仿真数据生成脚本
% 环境：20m x 40m 平地空间 (3D 算法仿真，高度常置为 0)
% 基站：默认 2 个，位于对角线 (0,0,0) 和 (20,40,0)
% 车辆：默认 4 辆，不同起点，运行 300 秒，包含加减速及最多 1 次转弯
% 采样：IMU 频率 100Hz (dt = 0.01s)，UWB 频率 10Hz (dt = 0.1s)
% =========================================================================

clc; clear; close all;

%% 1. 全局参数设置与动态数量声明
dt_imu = 0.01;              % IMU 采样时间 100Hz (0.01s)
dt_uwb = 0.1;               % UWB 采样时间 10Hz (0.1s)
t_end = 300;                % 运行时间 300秒

% 动态调整的基站和车辆数量 (满足后续动态调整需求)
Vehicle_num = 4;            
Anchor_num = 4;             

% 基站 3D 位置 (x, y, z) - 沿 20m x 40m 对角线放置于平地
anchors_pool = [ 0,  0, 0;
                20,  0,  0;
                0 , 40, 0;
                20, 40, 0]; 
% 截取指定数量的基站
anchors = anchors_pool(1:Anchor_num, :);

% 存储文件名与路径 (保存于当前 Data 文件夹下)
save_dir = 'E:\SE3_MLKF\Data'; 
trajectories_mat_name = sprintf('Trj_data_Veh%d_Anc%d_3D.mat', Vehicle_num, Anchor_num);

% 初始化总体轨迹结构体
trajectories = struct();
IMU_noise_params = struct();
UWB_noise_params = struct();

%% 2. 初始化车辆状态并生成轨迹 (包含加减速与至多 1 次转弯)
% 初始状态列表 [x, y, z, 初始速度, 初始航向角 theta]
init_states = [
    2,  5,  0, 0.05, pi/2;   % V1: 左下角，朝北
    16, 3,  0, 0.03, pi/2;   % V2: 右下角，朝北
    4,  35, 0, 0.05, 0;      % V3: 左上角，朝东
    18, 38, 0, 0.05, pi      % V4: 右上角，朝西
];

% 预设前4辆车的运动段配置 [加速度 a, 角速度 w, 持续时间 duration]
segs = cell(Vehicle_num, 1);
segs{1} = [
    0.001,  0,      60;      % 0-60s: 加速北行
    0,      0,      90;      % 60-150s: 匀速北行
    0,     -pi/60,  30;      % 150-180s: 右转90度向东 (30秒完成)
    0,      0,      60;      % 180-240s: 匀速东行
   -0.001,  0,      60       % 240-300s: 减速东行
];
segs{2} = [
    0.0012, 0,      50;      % 0-50s: 加速北行
    0,      0,      80;      % 50-130s: 匀速北行
    0,      pi/60,  30;      % 130-160s: 左转90度向西
    0,      0,      90;      % 160-250s: 匀速西行
   -0.0012, 0,      50       % 250-300s: 减速西行
];
segs{3} = [
    0.002,  0,      40;      % 0-40s: 加速东行
    0,      0,      60;      % 40-100s: 匀速东行
    0,     -pi/60,  30;      % 100-130s: 右转90度向南
    0,      0,      110;     % 130-240s: 匀速南行
   -0.001,  0,      60       % 240-300s: 减速南行
];
segs{4} = [
    0.001,  0,      40;      % 0-40s: 加速西行
    0,      0,      70;      % 40-110s: 匀速西行
    0,      pi/60,  30;      % 110-140s: 左转90度向南
    0,      0,      110;     % 140-250s: 匀速南行
   -0.001,  0,      50       % 250-300s: 减速南行
];

% 兼容车辆数变动时的自动拓展机制
for i = 5:Vehicle_num
    % 随机初始状态
    init_states(i, :) = [rand()*15+2, rand()*35+2, 0, 0.05, rand()*2*pi];
    segs{i} = [
        0.001,  0,      100;
        0,      0,      100;
       -0.001,  0,      100
    ];
end

% 循环生成每辆小车 100Hz 的连续真实轨迹
for i = 1:Vehicle_num
    v_name = sprintf('V%d', i);
    s0 = init_states(i, :);
    % 初始化
    trajectories.(v_name) = init_vehicle_3d(s0(1), s0(2), s0(3), s0(4), s0(5));
    
    % 拼接运行段
    v_seg = segs{i};
    for s = 1:size(v_seg, 1)
        trajectories.(v_name) = add_trajectory_segment_3d(...
            trajectories.(v_name), v_seg(s,1), v_seg(s,2), v_seg(s,3), dt_imu);
    end
    
    % 3D SO(3) 旋转矩阵序列 R_true 生成 (R 属于 SO(3))
    N_steps = length(trajectories.(v_name).Time_true);
    R_true = zeros(3, 3, N_steps);
    theta_arr = trajectories.(v_name).Theta_true;
    for k = 1:N_steps
        th = theta_arr(k);
        R_true(:, :, k) = [
            cos(th), -sin(th), 0;
            sin(th),  cos(th), 0;
            0,        0,       1
        ];
    end
    trajectories.(v_name).R_true = R_true;
end

%% 3. 生成 100Hz 的 3D IMU 数据 (包含 3D 偏置及重力影响)
% 设定 MEMS IMU 典型噪声与随机游走偏置参数
IMU_noise_params.sigma_na = 0.05;      % 加速度计白噪声标准差 (m/s^2)
IMU_noise_params.sigma_nw = 0.005;     % 陀螺仪白噪声标准差 (rad/s)
IMU_noise_params.sigma_ba = 0.002;     % 加速度偏置随机游走标准差 (m/s^2 * sqrt(s))
IMU_noise_params.sigma_bw = 0.0002;    % 陀螺仪偏置随机游走标准差 (rad/s * sqrt(s))

g_vec = [0; 0; -9.81];                 % 3D 导航系重力矢量

for i = 1:Vehicle_num
    v_name = sprintf('V%d', i);
    veh = trajectories.(v_name);
    t_imu = veh.Time_true;
    N_imu = length(t_imu);
    
    % 计算理想运动学下的 3D 世界系加速度
    ax_true = gradient(veh.Vx_true, dt_imu);
    ay_true = gradient(veh.Vy_true, dt_imu);
    az_true = gradient(veh.Vz_true, dt_imu); % 平地行驶高度不变，此项为 0
    
    % 将特定力转换到机体系：f_body = R^T * (a_world - g) [1, 5]
    a_body_ideal = zeros(N_imu, 3);
    for k = 1:N_imu
        R_k = veh.R_true(:, :, k);
        a_world = [ax_true(k); ay_true(k); az_true(k)];
        a_body_ideal(k, :) = (R_k' * (a_world - g_vec))';
    end
    
    % 计算理想角速度
    theta_unwrapped = unwrap(veh.Theta_true);
    wz_true = gradient(theta_unwrapped, dt_imu);
    omega_body_ideal = [zeros(N_imu, 2), wz_true]; % 仅绕 Z 轴旋转 [5]
    
    % 生成 3D 随机游走零偏 (Bias Evolution)
    ba_init = (rand(1, 3) - 0.5) * 0.1;    % 初始零偏
    bw_init = (rand(1, 3) - 0.5) * 0.01;
    b_a = ba_init + cumsum(randn(N_imu, 3) * IMU_noise_params.sigma_ba * sqrt(dt_imu), 1);
    b_w = bw_init + cumsum(randn(N_imu, 3) * IMU_noise_params.sigma_bw * sqrt(dt_imu), 1);
    
    % 引入高斯白噪声
    acc_noise = randn(N_imu, 3) * IMU_noise_params.sigma_na;
    gyro_noise = randn(N_imu, 3) * IMU_noise_params.sigma_nw;
    
    % 生成 IMU 测量值 (理想值 + 零偏 + 噪声)
    acc_m = a_body_ideal + b_a + acc_noise;
    gyro_m = omega_body_ideal + b_w + gyro_noise;
    
    % 写入结构体
    trajectories.(v_name).IMU_Time = t_imu;
    trajectories.(v_name).IMU_acc_m = acc_m;
    trajectories.(v_name).IMU_gyro_m = gyro_m;
    trajectories.(v_name).IMU_bias_a_true = b_a;
    trajectories.(v_name).IMU_bias_w_true = b_w;
end

%% 4. 生成 10Hz 的 UWB 3D 测距数据
UWB_noise_params.sigma_anc = 0.05;     % 基站测距白噪声标准差 (5cm 级别)
UWB_noise_params.sigma_rel = 0.07;     % 相对测距白噪声标准差 (7cm 级别)

% UWB 每10个 IMU 周期采样一次 (100Hz -> 10Hz 严格同步对齐)
idx_uwb = 1:10:N_imu;
t_uwb = t_imu(idx_uwb);
N_uwb = length(t_uwb);

% 提取 10Hz 时刻下所有车辆的 3D 真实位置
pos_true_uwb = zeros(N_uwb, 3, Vehicle_num);
for i = 1:Vehicle_num
    v_name = sprintf('V%d', i);
    pos_true_uwb(:, 1, i) = trajectories.(v_name).X_true(idx_uwb);
    pos_true_uwb(:, 2, i) = trajectories.(v_name).Y_true(idx_uwb);
    pos_true_uwb(:, 3, i) = trajectories.(v_name).Z_true(idx_uwb);
end

for i = 1:Vehicle_num
    v_name = sprintf('V%d', i);
    
    % ---------------- A. UWB 基站测距 ----------------
    % 数据格式: [时间, 基站1距离, 基站2距离, ...]
    UWB_Anchor = zeros(N_uwb, 1 + Anchor_num);
    UWB_Anchor(:, 1) = t_uwb;
    
    for a_idx = 1:Anchor_num
        dx = pos_true_uwb(:, 1, i) - anchors(a_idx, 1);
        dy = pos_true_uwb(:, 2, i) - anchors(a_idx, 2);
        dz = pos_true_uwb(:, 3, i) - anchors(a_idx, 3);
        dist_true = sqrt(dx.^2 + dy.^2 + dz.^2);
        
        % 叠加高斯噪声
        UWB_Anchor(:, 1 + a_idx) = dist_true + randn(N_uwb, 1) * UWB_noise_params.sigma_anc;
    end
    trajectories.(v_name).UWB_Anchor = UWB_Anchor;
    
    % ---------------- B. 车间相对测距 ----------------
    % 数据格式: [时间, V1距离, V2距离, V3距离, V4距离, ...]
    UWB_Relative = zeros(N_uwb, 1 + Vehicle_num);
    UWB_Relative(:, 1) = t_uwb;
    
    for j = 1:Vehicle_num
        if i == j
            UWB_Relative(:, 1 + j) = NaN; % 自身相对距离填充 NaN
        else
            dx = pos_true_uwb(:, 1, i) - pos_true_uwb(:, 1, j);
            dy = pos_true_uwb(:, 2, i) - pos_true_uwb(:, 2, j);
            dz = pos_true_uwb(:, 3, i) - pos_true_uwb(:, 3, j);
            dist_true = sqrt(dx.^2 + dy.^2 + dz.^2);
            
            UWB_Relative(:, 1 + j) = dist_true + randn(N_uwb, 1) * UWB_noise_params.sigma_rel;
        end
    end
    trajectories.(v_name).UWB_Relative = UWB_Relative;
end

%% 5. 结果可视化
figure('Name', 'Multi-Agent UWB/IMU 3D-Flat Trajectories', 'Position', [100, 100, 600, 800]);
hold on; grid on; axis equal;
xlim([-2, 22]); ylim([-2, 42]);
xlabel('X Position (m)'); ylabel('Y Position (m)');
title(sprintf('%d Vehicles Trajectories (300s, 20x40m Flat Space)', Vehicle_num));

% 绘制仿真边界
rectangle('Position', [0, 0, 20, 40], 'EdgeColor', 'k', 'LineWidth', 1.5, 'LineStyle', '--');

% 绘制基站位置
h_anchor = plot(anchors(:,1), anchors(:,2), '^', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'UWB Anchors');

colors = lines(Vehicle_num);
h_traj = zeros(1, Vehicle_num);

for i = 1:Vehicle_num
    v_name = sprintf('V%d', i);
    v_data = trajectories.(v_name);
   
    h_traj(i) = plot(v_data.X_true, v_data.Y_true, 'Color', colors(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Vehicle %d', i));

    % 绘制起始点和终点
    plot(v_data.X_true(1), v_data.Y_true(1), 'o', 'MarkerSize', 6, 'MarkerFaceColor', colors(i,:), 'Color', colors(i,:));
    plot(v_data.X_true(end), v_data.Y_true(end), '*', 'MarkerSize', 8, 'Color', colors(i,:));
end

legend([h_anchor, h_traj], 'Location', 'northeastoutside');
hold off;

%% 6. 数据导出
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
save_path = fullfile(save_dir, trajectories_mat_name);
save(save_path, 'trajectories', 'anchors', 'IMU_noise_params', 'UWB_noise_params', 'Vehicle_num', 'Anchor_num');
fprintf('仿真数据生成成功，已保存至：%s\n', save_path);

%% 辅助函数组

% 3D 小车状态初始化
function veh = init_vehicle_3d(x0, y0, z0, v0, theta0)
    veh.X_true = x0;
    veh.Y_true = y0;
    veh.Z_true = z0;
    veh.Vx_true = v0 * cos(theta0);
    veh.Vy_true = v0 * sin(theta0);
    veh.Vz_true = 0;
    veh.Theta_true = theta0;
    veh.Time_true = 0;
end

% 3D 理想运动状态迭代拼接
function veh = add_trajectory_segment_3d(veh, a, w, duration, dt)
    num_steps = round(duration / dt);
    
    curr_x = veh.X_true(end);
    curr_y = veh.Y_true(end);
    curr_z = veh.Z_true(end);
    
    curr_vx = veh.Vx_true(end);
    curr_vy = veh.Vy_true(end);
    curr_vz = veh.Vz_true(end);
    
    curr_theta = veh.Theta_true(end);
    curr_time = veh.Time_true(end);
    
    curr_v = sqrt(curr_vx^2 + curr_vy^2); % 当前地平面切向速标量
    
    new_X = zeros(num_steps, 1);
    new_Y = zeros(num_steps, 1);
    new_Z = zeros(num_steps, 1);
    new_Vx = zeros(num_steps, 1);
    new_Vy = zeros(num_steps, 1);
    new_Vz = zeros(num_steps, 1);
    new_Theta = zeros(num_steps, 1);
    new_Time = zeros(num_steps, 1);
    
    for i = 1:num_steps
        curr_time = curr_time + dt;
        
        % 位置积分
        curr_x = curr_x + curr_vx * dt;
        curr_y = curr_y + curr_vy * dt;
        curr_z = curr_z + curr_vz * dt;
        
        % 朝向积分
        curr_theta = curr_theta + w * dt;
        curr_theta = mod(curr_theta + pi, 2*pi) - pi;
        
        % 速度计算
        curr_v = curr_v + a * dt;
        curr_vx = curr_v * cos(curr_theta);
        curr_vy = curr_v * sin(curr_theta);
        
        new_X(i) = curr_x;
        new_Y(i) = curr_y;
        new_Z(i) = curr_z;
        new_Vx(i) = curr_vx;
        new_Vy(i) = curr_vy;
        new_Vz(i) = curr_vz;
        new_Theta(i) = curr_theta;
        new_Time(i) = curr_time;
    end
    
    veh.X_true = [veh.X_true; new_X];
    veh.Y_true = [veh.Y_true; new_Y];
    veh.Z_true = [veh.Z_true; new_Z];
    veh.Vx_true = [veh.Vx_true; new_Vx];
    veh.Vy_true = [veh.Vy_true; new_Vy];
    veh.Vz_true = [veh.Vz_true; new_Vz];
    veh.Theta_true = [veh.Theta_true; new_Theta];
    veh.Time_true = [veh.Time_true; new_Time];
end

