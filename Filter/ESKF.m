classdef ESKF < handle
    properties
        Vehicle_num         % 车辆数量 I
        states              % 结构体数组：.p (3x1), .v (3x1), .R (3x3) (加速度和角速度退化为输入)
        P                   % 联合协方差矩阵 [9I x 9I]
        Q_joint             % 联合过程噪声矩阵 [9I x 9I] (由IMU白噪声驱动)
        g_vec               % 3D重力加速度常数 [1]
    end
    
    methods
        function obj = ESKF(init_states_struct, init_P, IMU_noise_params)
            % ESKF 构造函数 (符合导师要求：将 IMU 视作输入，状态缩减至 9 维)
            %   init_states_struct: 初始状态结构体 (.p, .v, .R)
            %   init_P: 初始联合协方差 [9I x 9I]
            %   IMU_noise_params: 包含输入噪声 sigma_na, sigma_nw 的结构体
            
            obj.Vehicle_num = length(init_states_struct);
            obj.states = init_states_struct;
            obj.P = init_P;
            obj.g_vec = [0; 0; -9.81]; % 3D 重力矢量 [1]
            
            % 构建由输入白噪声驱动的离散系统过程噪声 Q_joint (9I x 9I)
            I = obj.Vehicle_num;
            obj.Q_joint = zeros(9*I, 9*I);
            
            % 经典输入噪声分配：位置块为0，速度块由加速度噪声驱动，姿态块由陀螺仪噪声驱动
            Q_single = diag([ ...
                1e-10 * ones(1, 3), ...                % 位置块过程噪声极小
                IMU_noise_params.sigma_na^2 * ones(1, 3), ... % 速度块过程噪声
                IMU_noise_params.sigma_nw^2 * ones(1, 3)  ... % 姿态块过程噪声
            ]);
            
            for idx = 1:I
                row_idx = (idx-1)*9 + (1:9);
                obj.Q_joint(row_idx, row_idx) = Q_single;
            end
        end
        
        function propagate(obj, imu_acc, imu_gyro, dt)
            % 标称状态与协方差时间传播 (符合导师要求：IMU 数据作为系统输入)
            % imu_acc:  [I x 3] 预校正后的加速度计输入 \tilde{a}_t [1]
            % imu_gyro: [I x 3] 预校正后的陀螺仪输入 \tilde{\omega}_t [1]
            % dt: 采样周期 \tau [1]
            
            I = obj.Vehicle_num;
            A_joint = zeros(9*I, 9*I);
            
            for i = 1:I
                p_t = obj.states(i).p;
                v_t = obj.states(i).v;
                R_t = obj.states(i).R;
                
                acc_tilde = imu_acc(i, :)';
                gyro_tilde = imu_gyro(i, :)';
                
                % 1. 经典名义状态积分传播 (IMU作为输入驱动位置、速度、姿态前向更新)
                acc_nav = R_t * acc_tilde + obj.g_vec; % 转换至导航系并叠加重力 [1]
                p_next = p_t + dt * v_t + 0.5 * dt^2 * acc_nav;
                v_next = v_t + dt * acc_nav;
                R_next = R_t * so3_exp(dt * gyro_tilde);
                
                obj.states(i).p = p_next;
                obj.states(i).v = v_next;
                obj.states(i).R = R_next;
                
                % 2. 构造 9 维误差状态转移矩阵 A_i (对应经典惯导误差方程)
                A_i = eye(9);
                A_i(1:3, 4:6) = dt * eye(3);                      % \delta p <- \delta v
                A_i(4:6, 7:9) = -dt * R_t * skew(acc_tilde);      % \delta v <- \delta \phi (比力驱动项)
                A_i(7:9, 7:9) = so3_exp(-dt * gyro_tilde);        % \delta \phi <- \delta \phi
                
                row_idx = (i-1)*9 + (1:9);
                A_joint(row_idx, row_idx) = A_i;
            end
            
            % 3. 联合协方差传播 (严格包含离散时间噪声 Q) [4]
            obj.P = A_joint * obj.P * A_joint' + obj.Q_joint * dt;
            obj.P = 0.5 * (obj.P + obj.P'); % 数值正定保护
        end
        
        function update(obj, anchors, uwb_anc, uwb_rel, sig_s, sig_z)
            % 经典紧耦合 ESKF 观测更新 (仅利用 UWB 测距修正惯导累积误差，IMU 不参与更新)
            % 输入:
            %   anchors: [K x 3] 基站位置 [1]
            %   uwb_anc: [I x K] 基站测距测量 (NaN表示无) [1]
            %   uwb_rel: [I x I] 车间相对测距测量 (NaN表示无) [1]
            
            I = obj.Vehicle_num;
            
            % --- 1. 活跃 UWB 测距统计 ---
            active_anc = [];   
            active_rel = [];   
            
            for i = 1:I
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
            
            M_anc = size(active_anc, 1);
            M_rel = size(active_rel, 1);
            M = M_anc + M_rel;
            
            if M == 0, return; end
            
            % 构造测量向量与 R 矩阵
            y_meas = [];
            R_list = [];
            for r = 1:M_anc
                y_meas = [y_meas; active_anc(r, 3)];
                R_list = [R_list; sig_s^2];
            end
            for r = 1:M_rel
                y_meas = [y_meas; active_rel(r, 3)];
                R_list = [R_list; sig_z^2];
            end
            R_cov = diag(R_list);
            
            % --- 2. 在当前标称预测位置计算残差与 9I 维观测雅可比 [4] ---
            h_val = zeros(M, 1);
            H_jac = zeros(M, 9*I); % 联合误差状态雅可比
            
            row_ptr = 1;
            
            % A. 基站测距雅可比 (仅关于位置误差 \delta p) [5]
            for r = 1:M_anc
                i = active_anc(r, 1); k = active_anc(r, 2);
                p_i = obj.states(i).p;
                c_k = anchors(k, :)';
                
                dist = norm(p_i - c_k);
                if dist < 1e-6, dist = 1e-6; end
                h_val(row_ptr) = dist;
                
                u_dir = (p_i - c_k)' / dist;
                H_jac(row_ptr, (i-1)*9 + (1:3)) = u_dir; % wrt \delta p_i
                row_ptr = row_ptr + 1;
            end
            
            % B. 车间测距雅可比 (关于 \delta p_i 和 \delta p_j) [5]
            for r = 1:M_rel
                i = active_rel(r, 1); j = active_rel(r, 2);
                p_i = obj.states(i).p;
                p_j = obj.states(j).p;
                
                dist = norm(p_i - p_j);
                if dist < 1e-6, dist = 1e-6; end
                h_val(row_ptr) = dist;
                
                u_dir = (p_i - p_j)' / dist;
                H_jac(row_ptr, (i-1)*9 + (1:3)) =  u_dir; % wrt \delta p_i
                H_jac(row_ptr, (j-1)*9 + (1:3)) = -u_dir; % wrt \delta p_j
                row_ptr = row_ptr + 1;
            end
            
            y_err = y_meas - h_val;
            
            % --- 3. 标准卡尔曼形式更新 ---
            S = H_jac * obj.P * H_jac' + R_cov;
            K = (obj.P * H_jac') / S;
            
            % 求解 9I 维误差状态更新向量 dx
            dx = K * y_err;
            
            % 更新状态协方差 P (对称性保护)
            obj.P = (eye(9*I) - K * H_jac) * obj.P;
            obj.P = 0.5 * (obj.P + obj.P');
            
            % --- 4. 标称状态 9维 流形回馈纠正 [7] ---
            for i = 1:I
                dx_i = dx((i-1)*9 + (1:9));
                
                % 对 9 维经典状态进行直和反馈修正
                obj.states(i).p = obj.states(i).p + dx_i(1:3);
                obj.states(i).v = obj.states(i).v + dx_i(4:6);
                obj.states(i).R = obj.states(i).R * so3_exp(dx_i(7:9));
            end
        end
    end
end