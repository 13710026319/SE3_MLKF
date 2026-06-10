classdef CMLKF < handle
    properties
        Vehicle_num         % 车辆数量 I
        states              % 结构体数组：.p (3x1), .v (3x1), .a (3x1), .R (3x3), .omega (3x1) [1]
        P                   % 联合协方差矩阵 [15I x 15I] [3]
        Q_joint             % 联合过程噪声矩阵 [15I x 15I] [4]
        g_vec               % 3D重力加速度常数 [1]
    end
    
    methods
        function obj = CMLKF(init_states_struct, init_P, Q_sigmas)
            % CMLKF 构造函数
            % Q_sigmas 包含以下过程噪声标准差：
            %   .sig_wp: 位置过程噪声标准差 (m) [1]
            %   .sig_wv: 速度过程噪声标准差 (m/s) [1]
            %   .sig_wa: 加速度过程噪声标准差 (m/s^2) [1]
            %   .sig_wR: 姿态过程噪声标准差 (rad) [1]
            %   .sig_womega: 角速度过程噪声标准差 (rad/s) [1]
            
            obj.Vehicle_num = length(init_states_struct);
            obj.states = init_states_struct;
            obj.P = init_P;
            obj.g_vec = [0; 0; -9.81]; % 导航系重力加速度 [1]
            
            % 构建联合过程噪声协方差矩阵 Q (15I x 15I)
            I = obj.Vehicle_num;
            obj.Q_joint = zeros(15*I, 15*I);
            
            Q_single = diag([ ...
                Q_sigmas.sig_wp^2 * ones(1, 3), ...
                Q_sigmas.sig_wv^2 * ones(1, 3), ...
                Q_sigmas.sig_wa^2 * ones(1, 3), ...
                Q_sigmas.sig_wR^2 * ones(1, 3), ...
                Q_sigmas.sig_womega^2 * ones(1, 3) ...
            ]);
            
            for idx = 1:I
                row_idx = (idx-1)*15 + (1:15);
                obj.Q_joint(row_idx, row_idx) = Q_single;
            end
        end
        
        function propagate(obj, dt)
            % 3.1 状态与协方差传播前向传播 (预测步，时间步长 dt = \tau) [3]
            % 依靠离散动力学进行预测，此架构不使用 IMU 作为前向输入的强约束
            
            I = obj.Vehicle_num;
            A_joint = zeros(15*I, 15*I);
            
            for i = 1:I
                % 1. 获取上一时刻 nominal 状态 [3]
                p_t = obj.states(i).p;
                v_t = obj.states(i).v;
                a_t = obj.states(i).a;
                R_t = obj.states(i).R;
                omega_t = obj.states(i).omega;
                
                % 2. 状态时间更新 (文档 Eq 17-21) [3]
                p_next = p_t + dt * v_t + 0.5 * dt^2 * a_t;
                v_next = v_t + dt * a_t;
                a_next = a_t;
                R_next = R_t * so3_exp(dt * omega_t); % 姿态矩阵传播 [3]
                omega_next = omega_t;
                
                % 保存预测值
                obj.states(i).p = p_next;
                obj.states(i).v = v_next;
                obj.states(i).a = a_next;
                obj.states(i).R = R_next;
                obj.states(i).omega = omega_next;
                
                % 3. 构造各车误差传播雅可比 A_i (文档 Eq 31) [4]
                I3 = eye(3);
                A_i = zeros(15, 15);
                A_i(1:3, 1:3)   = I3;
                A_i(1:3, 4:6)   = dt * I3;
                A_i(1:3, 7:9)   = 0.5 * dt^2 * I3;
                A_i(4:6, 4:6)   = I3;
                A_i(4:6, 7:9)   = dt * I3;
                A_i(7:9, 7:9)   = I3;
                
                % 姿态与角速度误差传播部分
                A_i(10:12, 10:12) = so3_exp(-dt * omega_t);        % exp(-\tau * [omega_t]x) [4]
                A_i(10:12, 13:15) = dt * so3_right_jacobian(dt * omega_t); % \tau * Jr(\tau * omega_t) [4]
                A_i(13:15, 13:15) = I3;
                
                row_idx = (i-1)*15 + (1:15);
                A_joint(row_idx, row_idx) = A_i;
            end
            
            % 4. 联合协方差更新 (文档 Eq 32) [4]
            obj.P = A_joint * obj.P * A_joint' + obj.Q_joint;
            obj.P = 0.5 * (obj.P + obj.P'); % 数值对称性保护
        end
        

        function update_imu_only(obj, imu_acc, imu_gyro, sig_acc, sig_gyro)
            % High-Frequency IMU-Only Update (新文档 3.2.4 节新增 EKF 更新)
            % 用于高频非UWB时刻，仅利用加速度计和陀螺仪测量进行状态修正 [1]
            
            I = obj.Vehicle_num;
            
            % 1. 构造联合测量残差 r_t 与 联合噪声协方差 R_IMU
            r_t = zeros(6*I, 1);
            H_IMU = zeros(6*I, 15*I);
            R_IMU = zeros(6*I, 6*I);
            
            for i = 1:I
                R_i = obj.states(i).R;
                a_i = obj.states(i).a;
                omega_i = obj.states(i).omega;
                
                acc_tilde = imu_acc(i, :)';
                gyro_tilde = imu_gyro(i, :)';
                
                % 计算单车 IMU 观测残差 (新文档 1.2 节测量方程) [1]
                r_a = acc_tilde - R_i' * (a_i - obj.g_vec);
                r_omega = gyro_tilde - omega_i;
                r_t((i-1)*6 + (1:6)) = [r_a; r_omega];
                
                % 2. 构造单车误差状态雅可比 H_IMU_t^i (6x15矩阵) [5]
                % 误差排布: [\delta p (1:3); \delta v (4:6); \delta a (7:9); \delta \phi (10:12); \delta \omega (13:15)] [2]
                H_i = zeros(6, 15);
                % 对应加速度测量方程对各误差状态求导 [5]
                H_i(1:3, 7:9)   = R_i';                                  % wrt \delta a
                H_i(1:3, 10:12) = skew(R_i' * (a_i - obj.g_vec));        % wrt \delta \phi
                % 对应角速度测量方程对各误差状态求导 [5]
                H_i(4:6, 13:15) = eye(3);                                % wrt \delta \omega
                
                H_IMU((i-1)*6 + (1:6), (i-1)*15 + (1:15)) = H_i;
                
                % 噪声协方差 R_IMU
                R_IMU((i-1)*6 + (1:6), (i-1)*6 + (1:6)) = diag([sig_acc^2 * ones(3, 1); sig_gyro^2 * ones(3, 1)]);
            end
            
            % 3. 标准卡尔曼形式更新 (无需非线性迭代优化，通过先验协方差矩阵 obj.P 进行正则化保护) [7]
            S = H_IMU * obj.P * H_IMU' + R_IMU;
            K = (obj.P * H_IMU') / S;
            
            % 计算 15I 维状态纠正偏差 Delta_theta
            dtheta = K * r_t;
            
            % 更新系统状态协方差 P (对称性保护)
            obj.P = (eye(15*I) - K * H_IMU) * obj.P;
            obj.P = 0.5 * (obj.P + obj.P');
            
            % 4. 标称状态 15维 流形修正 [7]
            for i = 1:I
                dx_i = dtheta((i-1)*15 + (1:15));
                obj.states(i) = manifold_add(obj.states(i), dx_i);
            end
        end


        function update(obj, imu_acc, imu_gyro, anchors, uwb_anc, uwb_rel, ...
                        sig_acc, sig_gyro, sig_s, sig_z)
            % 3.2 & 3.3 集中式观测更新步骤 [4, 6]
            % 输入:
            %   imu_acc, imu_gyro: [I x 3] 去偏置后的 IMU 测量值 (已扣除 biases) [1]
            %   anchors:           [K x 3] UWB基站位置 [1]
            %   uwb_anc:           [I x K] 基站测距测量值 (NaN表示无) [1]
            %   uwb_rel:           [I x I] 车辆间相对测距测量值 (NaN表示无) [1]
            
            I = obj.Vehicle_num;
            
            % --- 1. 活跃观测链路统计 ---
            active_acc = [];   % [vehicle_idx, val_x, val_y, val_z]
            active_gyro = [];  % [vehicle_idx, val_x, val_y, val_z]
            active_anc = [];   % [vehicle_idx, anchor_idx, val]
            active_rel = [];   % [veh_i, veh_j, val]
            
            for i = 1:I
                if ~any(isnan(imu_acc(i, :)))
                    active_acc = [active_acc; i, imu_acc(i, :)];
                end
                if ~any(isnan(imu_gyro(i, :)))
                    active_gyro = [active_gyro; i, imu_gyro(i, :)];
                end
                for k = 1:size(anchors, 1)
                    if ~isnan(uwb_anc(i, k))
                        active_anc = [active_anc; i, k, uwb_anc(i, k)];
                    end
                end
                for j = (i+1):I
                    if ~isnan(uwb_rel(i, j))
                        active_rel = [active_rel; i, j, uwb_rel(i, j)];
                    end
                end
            end
            
            M_acc  = size(active_acc, 1);
            M_gyro = size(active_gyro, 1);
            M_anc  = size(active_anc, 1);
            M_rel  = size(active_rel, 1);
            M = 3*M_acc + 3*M_gyro + M_anc + M_rel; % 总测量维度
            
            if M == 0, return; end
            
            % 构造测量向量 y_meas
            y_meas = [];
            R_list = [];
            for r = 1:M_acc
                y_meas = [y_meas; active_acc(r, 2:4)'];
                R_list = [R_list; sig_acc^2 * ones(3, 1)];
            end
            for r = 1:M_gyro
                y_meas = [y_meas; active_gyro(r, 2:4)'];
                R_list = [R_list; sig_gyro^2 * ones(3, 1)];
            end
            for r = 1:M_anc
                y_meas = [y_meas; active_anc(r, 3)];
                R_list = [R_list; sig_s^2];
            end
            for r = 1:M_rel
                y_meas = [y_meas; active_rel(r, 3)];
                R_list = [R_list; sig_z^2];
            end
            inv_R = diag(1 ./ R_list);
            
            % 保存先验 nominal 状态以备后用
            states_prior = obj.states;
            
            % --- 2. 迭代高斯-牛顿优化 (仅对测量流形 \chi 优化, 维度为 12I) [4, 5] ---
            % 提取 \chi 优化初值 (Eq 15) [2]
            chi_opt = struct('p', {}, 'a', {}, 'R', {}, 'omega', {});
            for i = 1:I
                chi_opt(i).p     = obj.states(i).p;
                chi_opt(i).a     = obj.states(i).a;
                chi_opt(i).R     = obj.states(i).R;
                chi_opt(i).omega = obj.states(i).omega;
            end
            
            max_iter = 10;
            tol = 1e-4;
            
            for iter = 1:max_iter
                % 评估当前估计下的 h(\chi) 与 雅可比 H
                h_val = zeros(M, 1);
                H_jac = zeros(M, 12*I);
                
                row_ptr = 1;
                
                % A. 加速度计测距 Jacobian [5]
                for r = 1:M_acc
                    i = active_acc(r, 1);
                    R_i = chi_opt(i).R;
                    a_i = chi_opt(i).a;
                    
                    % 预测测量值 (文档 Eq 7) [1]
                    h_val(row_ptr:row_ptr+2) = R_i' * (a_i - obj.g_vec);
                    
                    % 局部扰动导数 (文档 Eq 37) [5]
                    H_jac(row_ptr:row_ptr+2, (i-1)*12 + (4:6)) = R_i'; % wrt \delta a
                    H_jac(row_ptr:row_ptr+2, (i-1)*12 + (7:9)) = skew(R_i' * (a_i - obj.g_vec)); % wrt \delta \phi
                    row_ptr = row_ptr + 3;
                end
                
                % B. 陀螺仪测距 Jacobian [5]
                for r = 1:M_gyro
                    i = active_gyro(r, 1);
                    omega_i = chi_opt(i).omega;
                    
                    h_val(row_ptr:row_ptr+2) = omega_i;
                    
                    % 局部扰动导数 (文档 Eq 38) [5]
                    H_jac(row_ptr:row_ptr+2, (i-1)*12 + (10:12)) = eye(3); % wrt \delta \omega
                    row_ptr = row_ptr + 3;
                end
                
                % C. UWB 基站测距 Jacobian [5]
                for r = 1:M_anc
                    i = active_anc(r, 1);
                    k = active_anc(r, 2);
                    p_i = chi_opt(i).p;
                    c_k = anchors(k, :)';
                    
                    dist = norm(p_i - c_k);
                    if dist < 1e-6, dist = 1e-6; end
                    h_val(row_ptr) = dist;
                    
                    % 局部扰动导数 (文档 Eq 39) [5]
                    u_dir = (p_i - c_k)' / dist;
                    H_jac(row_ptr, (i-1)*12 + (1:3)) = u_dir; % wrt \delta p
                    row_ptr = row_ptr + 1;
                end
                
                % D. UWB 车间测距 Jacobian [5]
                for r = 1:M_rel
                    i = active_rel(r, 1);
                    j = active_rel(r, 2);
                    p_i = chi_opt(i).p;
                    p_j = chi_opt(j).p;
                    
                    dist = norm(p_i - p_j);
                    if dist < 1e-6, dist = 1e-6; end
                    h_val(row_ptr) = dist;
                    
                    % 局部扰动导数 (文档 Eq 40) [5]
                    u_dir = (p_i - p_j)' / dist;
                    H_jac(row_ptr, (i-1)*12 + (1:3)) = u_dir;   % wrt \delta p_i
                    H_jac(row_ptr, (j-1)*12 + (1:3)) = -u_dir;  % wrt \delta p_j
                    row_ptr = row_ptr + 1;
                end
                
                residual = y_meas - h_val;
                Hessian_pos = H_jac' * inv_R * H_jac + 1e-4 * eye(12*I); % Levenberg-Marquardt阻尼 [5]
                Delta_s = Hessian_pos \ (H_jac' * inv_R * residual);     % 解局部增量 Eq 41 [5]
                
                % 本地流形退回/更新 (由于排除了速度，用 12D 退回形式更新) [5]
                for i = 1:I
                    ds_i = Delta_s((i-1)*12 + (1:12));
                    chi_opt(i).p     = chi_opt(i).p + ds_i(1:3);
                    chi_opt(i).a     = chi_opt(i).a + ds_i(4:6);
                    chi_opt(i).R     = chi_opt(i).R * so3_exp(ds_i(7:9));
                    chi_opt(i).omega = chi_opt(i).omega + ds_i(10:12);
                end
                
                if norm(Delta_s) < tol
                    break;
                end
            end
            
            % GN 结束，在收敛最优值 chi_L 处重新严谨估计最后一版 H_jac [10]
            chi_L = chi_opt;
            H_jac_ML = zeros(M, 12*I);
            row_ptr = 1;
            for r = 1:M_acc
                i = active_acc(r, 1);
                H_jac_ML(row_ptr:row_ptr+2, (i-1)*12 + (4:6)) = chi_L(i).R';
                H_jac_ML(row_ptr:row_ptr+2, (i-1)*12 + (7:9)) = skew(chi_L(i).R' * (chi_L(i).a - obj.g_vec));
                row_ptr = row_ptr + 3;
            end
            for r = 1:M_gyro
                i = active_gyro(r, 1);
                H_jac_ML(row_ptr:row_ptr+2, (i-1)*12 + (10:12)) = eye(3);
                row_ptr = row_ptr + 3;
            end
            for r = 1:M_anc
                i = active_anc(r, 1); k = active_anc(r, 2);
                u_dir = (chi_L(i).p - anchors(k, :)')' / norm(chi_L(i).p - anchors(k, :)');
                H_jac_ML(row_ptr, (i-1)*12 + (1:3)) = u_dir;
                row_ptr = row_ptr + 1;
            end
            for r = 1:M_rel
                i = active_rel(r, 1); j = active_rel(r, 2);
                u_dir = (chi_L(i).p - chi_L(j).p)' / norm(chi_L(i).p - chi_L(j).p);
                H_jac_ML(row_ptr, (i-1)*12 + (1:3)) = u_dir;
                H_jac_ML(row_ptr, (j-1)*12 + (1:3)) = -u_dir;
                row_ptr = row_ptr + 1;
            end
            
            % 提取似然等效信息矩阵 (文档 Eq 44) [6]
            inv_Xi_t = H_jac_ML' * inv_R * H_jac_ML;
            
            % 统计提取与空间雅可比融合 (文档 3.2.4 & 3.3 节) --- [6]
            mu_t = zeros(12*I, 1);
            J_t = zeros(12*I, 12*I);
            
            for i = 1:I
                % 获取该车的先验预测值
                p_hat = states_prior(i).p;
                a_hat = states_prior(i).a;
                R_hat = states_prior(i).R;
                omega_hat = states_prior(i).omega;
                
                % A. 计算该车在测量流形上的 ML 状态偏差 mu_t (文档 Eq 43) [6]
                mu_t((i-1)*12 + (1:3))   = chi_L(i).p - p_hat;
                mu_t((i-1)*12 + (4:6))   = chi_L(i).a - a_hat;
                mu_t((i-1)*12 + (7:9))   = so3_log(R_hat' * chi_L(i).R);
                mu_t((i-1)*12 + (10:12)) = chi_L(i).omega - omega_hat;
                
                % B. 计算该车对应的逆右雅可比分量 Jr_inv (使用 -mu_t 代替 sc_0) [6, 7]
                phi_e_i = -mu_t((i-1)*12 + (7:9)); % 旋转偏差，严格等于 so3_log(chi_L(i).R' * R_hat)
                Jr_inv_i = so3_inv_right_jacobian(phi_e_i);
                
                % 构造单车空间雅可比变换块 J_t^i (文档 Eq 53) [7]
                J_i_t = eye(12);
                J_i_t(7:9, 7:9) = Jr_inv_i;
                
                row_idx = (i-1)*12 + (1:12);
                J_t(row_idx, row_idx) = J_i_t;
            end
            
            % 计算先验切平面的 mapped 测量信息矩阵 (文档 Eq 57) [7]
            inv_Xi_t_prior = J_t' * inv_Xi_t * J_t;
            
            % 计算先验切平面 mapped 信息向量 (已修正 2 倍过大 Bug) [7]
            eta_t_prior = J_t' * inv_Xi_t * mu_t;

            % --- 4. 误差切空间降维选择与融合更新 (文档 Eq 59-62) --- [7]
            % 构造切平面状态选择矩阵 \pi (12I x 15I) 排除速度 [2, 7]
            pi_mat = zeros(12*I, 15*I);
            for i = 1:I
                pi_i = zeros(12, 15);
                pi_i(1:3, 1:3)     = eye(3);  % 选择 \delta p
                pi_i(4:6, 7:9)     = eye(3);  % 选择 \delta a
                pi_i(7:9, 10:12)   = eye(3);  % 选择 \delta \phi
                pi_i(10:12, 13:15) = eye(3);  % 选择 \delta \omega
                
                row_idx = (i-1)*12 + (1:12);
                col_idx = (i-1)*15 + (1:15);
                pi_mat(row_idx, col_idx) = pi_i;
            end
            
            % 降维投影信息项 (文档 Eq 59-60) [7]
            Lambda_t = pi_mat' * inv_Xi_t_prior * pi_mat;
            lambda_t = pi_mat' * eta_t_prior;
            
            % 应用 Woodbury 矩阵求逆引理更新系统后验误差协方差，免去显式 inv 带来的奇异危险 [7]
            Sigma_prior = obj.P;
            S_eff = eye(15*I) + Lambda_t * Sigma_prior;
            Sigma_post = Sigma_prior - Sigma_prior * (S_eff \ (Lambda_t * Sigma_prior));
            Sigma_post = 0.5 * (Sigma_post + Sigma_post'); % 对称正定保护
            obj.P = Sigma_post;
            
            % 计算流形反撤更新增量 (文档 Eq 62) [7]
            dx_joint = Sigma_post * lambda_t;
            
            % 局部误差增量回馈叠加至名义状态 [7]
            for i = 1:I
                dx_i = dx_joint((i-1)*15 + (1:15));
                obj.states(i) = manifold_add(obj.states(i), dx_i);
            end
        end
    end
end
