% =========================================================================
% generate_data_3D.m (平滑物理连续版)
% 3D多车(无人机行为) UWB/IMU 协同定位仿真数据生成脚本
% 空间范围：20m x 40m x 8m (高度严格控制在 2m ~ 7m)
% 改进重点：为转弯设计 3秒 平滑原地过渡段，彻底消除瞬间转向带来的数值发散
% =========================================================================

clc; clear; close all;

%% 1. 全局参数设置与保存路径
dt_imu = 0.01;              % IMU 采样时间 100Hz (0.01s)
dt_uwb = 0.1;               % UWB 采样时间 10Hz (0.1s)
t_end = 300;                % 运行时间 300秒
N_steps = round(t_end / dt_imu) + 1; % 30001 个采样点

Vehicle_num = 4;            % 车辆数量
Anchor_num = 4;             % 基站数量

% 基站部署 (3D: x, y, z) - 位于四个角落，高度不相等，打破共面奇异
anchors = [ 0,  0, 1.5;     % 角落 1
           20,  0, 5.5;     % 角落 2
           20, 40, 3.0;     % 角落 3
            0, 40, 4.0];    % 角落 4

% 保存路径
save_dir = 'E:\SE3_MLKF\Data'; 
trajectories_mat_name = sprintf('Trj_data_Veh%d_Anc%d_3D.mat', Vehicle_num, Anchor_num);

% 噪声参数
IMU_noise_params.sigma_na = 0.05;      
IMU_noise_params.sigma_nw = 0.005;     
IMU_noise_params.sigma_ba = 0.002;     
IMU_noise_params.sigma_bw = 0.0002;    

UWB_noise_params.sigma_anc = 0.05;     
UWB_noise_params.sigma_rel = 0.07;     

%% 2. 独立车辆动力学规划 (含 3秒 平滑原地转弯)
init_configs = [
    5,  8,  3.0,  0.15,  pi/2;   % V1: 高度 3.0m，朝北
   16,  4,  5.0,  0.156,  pi/2;  % V2: 高度 5.0m，朝北
    2, 33,  6.0,  0.14,  0;      % V3: 高度 6.0m，朝东
   18, 36,  2.5,  0.15,  pi;     % V4: 高度 2.5m，朝西
];

trajectories = struct();

for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    cfg = init_configs(n, :);
    
    P_true = zeros(N_steps, 3);
    V_true = zeros(N_steps, 3);
    A_true = zeros(N_steps, 3);
    Theta_true = zeros(N_steps, 1);
    
    % 初始条件
    P_true(1, :) = cfg(1:3);
    v_mag = cfg(4);
    th_curr = cfg(5);
    V_true(1, :) = [v_mag * cos(th_curr), v_mag * sin(th_curr), 0];
    Theta_true(1) = th_curr;
    
    for k = 2:N_steps
        t = (k-1) * dt_imu;
        
        a_curr = [0; 0; 0];
        th_curr = Theta_true(k-1);
        
        switch n
            case 1 % Vehicle 1 (朝北 -> 减速 -> 平滑右转90度 -> 爬升 -> 加速东行)
                if t <= 60
                    a_curr = [0; 0; 0];
                elseif t > 60 && t <= 70 % 减速至静止
                    a_curr = [0; -0.015; 0];
                elseif t > 70 && t <= 73 % 原地平滑右转 90度 (北 -> 东)
                    a_curr = [0; 0; 0];
                    th_curr = pi/2 - (pi/2)/3 * (t - 70); % 【已修正】3s内均匀旋转至0
                elseif t > 73 && t <= 83 % 原地加速攀升 (3.0m -> 4.5m)
                    th_curr = 0;
                    a_curr = [0; 0; 0.015];
                elseif t > 83 && t <= 93 % 原地攀升减速至静止
                    th_curr = 0;
                    a_curr = [0; 0; -0.015];
                elseif t > 93 && t <= 102 % 东向加速1
                    th_curr = 0;
                    a_curr = [0.002; 0; 0];
                elseif t > 102 && t <= 117 % 东向加速2
                    th_curr = 0;
                    a_curr = [0.003; 0; 0];
                else % 匀速运行
                    th_curr = 0;
                    a_curr = [0; 0; 0];
                end
                
            case 2 % Vehicle 2 (朝北 -> 减速 -> 平滑左转90度 -> 下降 -> 加速西行)
                if t <= 50
                    a_curr = [0; 0; 0];
                elseif t > 50 && t <= 62 % 减速至静止
                    a_curr = [0; -0.013; 0];
                elseif t > 62 && t <= 65 % 原地平滑左转 90度 (北 -> 西)
                    a_curr = [0; 0; 0];
                    th_curr = pi/2 + (pi/2)/3 * (t - 62); % 【已修正】3s内均匀旋转至pi
                elseif t > 65 && t <= 74 % 原地加速下降 (5.0m -> 3.5m)
                    th_curr = pi;
                    a_curr = [0; 0; -0.0185];
                elseif t > 74 && t <= 83 % 原地下降减速至静止
                    th_curr = pi;
                    a_curr = [0; 0; 0.0185];
                elseif t > 83 && t <= 97 % 西向加速1
                    th_curr = pi;
                    a_curr = [-0.0013; 0; 0];
                elseif t > 97 && t <= 117 % 西向加速2
                    th_curr = pi;
                    a_curr = [-0.0025; 0; 0];
                else % 匀速运行
                    th_curr = pi;
                    a_curr = [0; 0; 0];
                end
                
            case 3 % Vehicle 3 (朝东 -> 减速 -> 平滑右转90度 -> 下降 -> 加速南行)
                if t <= 70
                    a_curr = [0; 0; 0];
                elseif t > 70 && t <= 80 % 减速至静止
                    a_curr = [-0.014; 0; 0];
                elseif t > 80 && t <= 83 % 原地平滑右转 90度 (东 -> 南)
                    a_curr = [0; 0; 0];
                    th_curr = 0 - (pi/2)/3 * (t - 80); % 【已修正】3s内均匀旋转至-pi/2
                elseif t > 83 && t <= 91 % 原地加速下降 (6.0m -> 4.0m)
                    th_curr = -pi/2;
                    a_curr = [0; 0; -0.03125];
                elseif t > 91 && t <= 99 % 原地下降减速至静止
                    th_curr = -pi/2;
                    a_curr = [0; 0; 0.03125];
                elseif t > 99 && t <= 112 % 南向加速1
                    th_curr = -pi/2;
                    a_curr = [0; -0.002; 0];
                elseif t > 112 && t <= 127 % 南向加速2
                    th_curr = -pi/2;
                    a_curr = [0; -0.003; 0];
                else % 匀速运行
                    th_curr = -pi/2;
                    a_curr = [0; 0; 0];
                end
                
            case 4 % Vehicle 4 (朝西 -> 减速 -> 平滑左转90度 -> 爬升 -> 加速南行)
                if t <= 60
                    a_curr = [0; 0; 0];
                elseif t > 60 && t <= 70 % 减速至静止
                    a_curr = [0.015; 0; 0];
                elseif t > 70 && t <= 73 % 原地平滑左转 90度 (西 -> 南)
                    a_curr = [0; 0; 0];
                    th_curr = pi + (pi/2)/3 * (t - 70); % 【已修正】3s内均匀旋转至3*pi/2
                elseif t > 73 && t <= 83 % 原地加速爬升 (2.5m -> 4.5m)
                    th_curr = 3*pi/2;
                    a_curr = [0; 0; 0.02];
                elseif t > 83 && t <= 93 % 原地爬升减速至静止
                    th_curr = 3*pi/2;
                    a_curr = [0; 0; -0.02];
                elseif t > 93 && t <= 102 % 南向加速1
                    th_curr = 3*pi/2;
                    a_curr = [0; -0.002; 0];
                elseif t > 102 && t <= 117 % 南向加速2
                    th_curr = 3*pi/2;
                    a_curr = [0; -0.004; 0];
                else % 匀速运行
                    th_curr = 3*pi/2;
                    a_curr = [0; 0; 0];
                end
        end
        
        % 欧拉积分更新速度与位置 (保证物理严格连续) [1]
        A_true(k-1, :) = a_curr';
        V_true(k, :)   = V_true(k-1, :) + A_true(k-1, :) * dt_imu;
        P_true(k, :)   = P_true(k-1, :) + V_true(k-1, :) * dt_imu + 0.5 * A_true(k-1, :) * dt_imu^2;
        Theta_true(k)  = th_curr;
    end
    A_true(end, :) = [0, 0, 0]; 
    
    % 保存真值分量 [1]
    trajectories.(v_name).Time_true = (0:dt_imu:t_end)';
    trajectories.(v_name).X_true = P_true(:, 1);
    trajectories.(v_name).Y_true = P_true(:, 2);
    trajectories.(v_name).Z_true = P_true(:, 3);
    trajectories.(v_name).Vx_true = V_true(:, 1);
    trajectories.(v_name).Vy_true = V_true(:, 2);
    trajectories.(v_name).Vz_true = V_true(:, 3);
    trajectories.(v_name).Theta_true = Theta_true;
    
    % 旋转矩阵序列 R_true [1, 2]
    R_true = zeros(3, 3, N_steps);
    for k = 1:N_steps
        th = Theta_true(k);
        R_true(:, :, k) = [
            cos(th), -sin(th), 0;
            sin(th),  cos(th), 0;
            0,        0,       1
        ];
    end
    trajectories.(v_name).R_true = R_true;
end

%% 3. 生成含有零偏与白噪声的 100Hz 3D IMU 测量信号
for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    veh = trajectories.(v_name);
    
    % 计算理想 3D 角速度 (绕 Z 轴旋转，此时因为2秒平滑过渡，wz极其平滑、合理)
    theta_unwrapped = unwrap(veh.Theta_true);
    wz_true = gradient(theta_unwrapped, dt_imu);
    omega_body_ideal = [zeros(N_steps, 2), wz_true];
    
    % 计算理想特定力：f = R^T * (a - g) [1, 5]
    a_body_ideal = zeros(N_steps, 3);
    g_vec = [0; 0; -9.81];
    for k = 1:N_steps
        R_k = veh.R_true(:, :, k);
        a_world = [A_true(k, 1); A_true(k, 2); A_true(k, 3)];
        a_body_ideal(k, :) = (R_k' * (a_world - g_vec))';
    end
    
    % 模拟高保真零偏 (Bias) 游走
    ba_init = (rand(1, 3) - 0.5) * 0.1;    
    bw_init = (rand(1, 3) - 0.5) * 0.01;
    b_a = ba_init + cumsum(randn(N_steps, 3) * IMU_noise_params.sigma_ba * sqrt(dt_imu), 1);
    b_w = bw_init + cumsum(randn(N_steps, 3) * IMU_noise_params.sigma_bw * sqrt(dt_imu), 1);
    
    % 叠加噪声 [1]
    acc_noise = randn(N_steps, 3) * IMU_noise_params.sigma_na;
    gyro_noise = randn(N_steps, 3) * IMU_noise_params.sigma_nw;
    
    acc_m = a_body_ideal + b_a + acc_noise;
    gyro_m = omega_body_ideal + b_w + gyro_noise;
    
    % 保存 IMU 数据
    trajectories.(v_name).IMU_Time = veh.Time_true;
    trajectories.(v_name).IMU_acc_m = acc_m;
    trajectories.(v_name).IMU_gyro_m = gyro_m;
    trajectories.(v_name).IMU_bias_a_true = b_a;
    trajectories.(v_name).IMU_bias_w_true = b_w;
end

%% 4. 生成 10Hz 同步 UWB 测距信号
idx_uwb = 1:10:N_steps;
t_uwb = trajectories.V1.Time_true(idx_uwb);
N_uwb = length(t_uwb);

pos_true_uwb = zeros(N_uwb, 3, Vehicle_num);
for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    pos_true_uwb(:, 1, n) = trajectories.(v_name).X_true(idx_uwb);
    pos_true_uwb(:, 2, n) = trajectories.(v_name).Y_true(idx_uwb);
    pos_true_uwb(:, 3, n) = trajectories.(v_name).Z_true(idx_uwb);
end

for n = 1:Vehicle_num
    v_name = sprintf('V%d', n);
    
    % A. UWB 基站 3D 测距 [1]
    UWB_Anchor = zeros(N_uwb, 1 + Anchor_num);
    UWB_Anchor(:, 1) = t_uwb;
    for a_idx = 1:Anchor_num
        dx = pos_true_uwb(:, 1, n) - anchors(a_idx, 1);
        dy = pos_true_uwb(:, 2, n) - anchors(a_idx, 2);
        dz = pos_true_uwb(:, 3, n) - anchors(a_idx, 3);
        dist_true = sqrt(dx.^2 + dy.^2 + dz.^2);
        
        UWB_Anchor(:, 1 + a_idx) = dist_true + randn(N_uwb, 1) * UWB_noise_params.sigma_anc;
    end
    trajectories.(v_name).UWB_Anchor = UWB_Anchor;
    
    % B. UWB 3D 车间相对测距 [1]
    UWB_Relative = zeros(N_uwb, 1 + Vehicle_num);
    UWB_Relative(:, 1) = t_uwb;
    for j = 1:Vehicle_num
        if n == j
            UWB_Relative(:, 1 + j) = NaN;
        else
            dx = pos_true_uwb(:, 1, n) - pos_true_uwb(:, 1, j);
            dy = pos_true_uwb(:, 2, n) - pos_true_uwb(:, 2, j);
            dz = pos_true_uwb(:, 3, n) - pos_true_uwb(:, 3, j);
            dist_true = sqrt(dx.^2 + dy.^2 + dz.^2);
            
            UWB_Relative(:, 1 + j) = dist_true + randn(N_uwb, 1) * UWB_noise_params.sigma_rel;
        end
    end
    trajectories.(v_name).UWB_Relative = UWB_Relative;
end

%% 5. 轨迹可视化 (三维立体作图以确认高度不共面)
figure('Name', 'Multi-Agent 3D Trajectories (Uncoplanar Environment)', 'Position', [100, 100, 800, 600]);
hold on; grid on; axis equal;
xlabel('X Position (m)'); ylabel('Y Position (m)'); zlabel('Height Z (m)');
title('4 Vehicles 3D Trajectories & Anchor Deployment (Altitude Envelope: 2m-7m)');

% 绘制边界立方体
line([0, 20, 20, 0, 0], [0, 0, 40, 40, 0], [0, 0, 0, 0, 0], 'Color', 'k', 'LineStyle', '--');
line([0, 20, 20, 0, 0], [0, 0, 40, 40, 0], [8, 8, 8, 8, 8], 'Color', 'k', 'LineStyle', '--');
for corner = [0, 20]
    for side = [0, 40]
        line([corner, corner], [side, side], [0, 8], 'Color', 'k', 'LineStyle', '--');
    end
end

% 绘制不共面基站
h_anchor = plot3(anchors(:,1), anchors(:,2), anchors(:,3), '^', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'DisplayName', 'Uncoplanar Anchors');

colors = lines(Vehicle_num);
h_traj = zeros(1, Vehicle_num);

for n = 1:Vehicle_num
    v_data = trajectories.(sprintf('V%d', n));
    h_traj(n) = plot3(v_data.X_true, v_data.Y_true, v_data.Z_true, 'Color', colors(n,:), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Vehicle %d', n));
    
    % 标记起始点和终点
    plot3(v_data.X_true(1), v_data.Y_true(1), v_data.Z_true(1), 'o', 'MarkerSize', 8, 'MarkerFaceColor', colors(n,:), 'Color', colors(n,:));
    plot3(v_data.X_true(end), v_data.Y_true(end), v_data.Z_true(end), '*', 'MarkerSize', 10, 'Color', colors(n,:));
end

view(3); 
legend([h_anchor, h_traj], 'Location', 'northeastoutside');
hold off;

%% 6. 数据保存
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
save_path = fullfile(save_dir, trajectories_mat_name);
save(save_path, 'trajectories', 'anchors', 'IMU_noise_params', 'UWB_noise_params', 'Vehicle_num', 'Anchor_num');
fprintf('3D 仿真数据生成成功，已成功保存至：%s\n', save_path);