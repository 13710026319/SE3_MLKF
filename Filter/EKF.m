classdef EKF < handle
    properties
        Vehicle_num         % 车辆数量 I
        states              % 结构体数组：.p (3x1), .v (3x1), .a (3x1), .R (3x3), .omega (3x1) [1]
        P                   % 联合协方差矩阵 [15I x 15I] [3]
        Q_joint             % 联合过程噪声矩阵 [15I x 15I] [4]
        g_vec               % 3D重力加速度常数 [1]
    end
    
    methods
        function obj = EKF(init_states_struct, init_P, Q_sigmas)
            % EKF 构造函数 (参数与 CMLKF 完全一致，确保公平对比)
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
            % 3.1 状态与协方差传播前向传播 (预测步与 CMLKF 保持完全一致) [3]
            I = obj.Vehicle_num;
            A_joint = zeros(15*I, 15*I);
            
            for i = 1:I
                p_t = obj.states(i).p;
                v_t = obj.states(i).v;
                a_t = obj.states(i).a;
                R_t = obj.states(i).R;
                omega_t = obj.states(i).omega;
                
                % 状态时间更新 (文档 Eq 17-21) [3]
                p_next = p_t + dt * v_t + 0.5 * dt^2 * a_t;
                v_next = v_t + dt * a_t;
                a_next = a_t;
                R_next = R_t * so3_exp(dt * omega_t);
                omega_next = omega_t;
                
                obj.states(i).p = p_next;
                obj.states(i).v = v_next;
                obj.states(i).a = a_next;
                obj.states(i).R = R_next;
                obj.states(i).omega = omega_next;
                
                % 构造各车误差传播雅可比 A_i (文档 Eq 31) [4]
                I3 = eye(3);
                A_i = zeros(15, 15);
                A_i(1:3, 1:3)   = I3;
                A_i(1:3, 4:6)   = dt * I3;
                A_i(1:3, 7:9)   = 0.5 * dt^2 * I3;
                A_i(4:6, 4:6)   = I3;
                A_i(4:6, 7:9)   = dt * I3;
                A_i(7:9, 7:9)   = I3;
                
                A_i(10:12, 10:12) = so3_exp(-dt * omega_t);
                A_i(10:12, 13:15) = dt * so3_right_jacobian(dt * omega_t);
                A_i(13:15, 13:15) = I3;
                
                row_idx = (i-1)*15 + (1:15);
                A_joint(row_idx, row_idx) = A_i;
            end
            
            % 联合协方差更新 (严格对齐公式 32) [4]
            obj.P = A_joint * obj.P * A_joint' + obj.Q_joint;
            obj.P = 0.5 * (obj.P + obj.P'); 
        end
        
        function update(obj, imu_acc, imu_gyro, anchors, uwb_anc, uwb_rel, ...
                        sig_acc, sig_gyro, sig_s, sig_z)
            % 经典集中式误差状态 EKF (ESEKF) 观测更新步
            % 直接基于先验预测状态进行单次线性化更新，无 GN 迭代 [4]
            
            I = obj.Vehicle_num;
            
            % --- 1. 活跃观测链路统计 (与 CMLKF 保持完全相同的结构) ---
            active_acc = [];   
            active_gyro = [];  
            active_anc = [];   
            active_rel = [];   
            
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
            M = 3*M_acc + 3*M_gyro + M_anc + M_rel; 
            
            if M == 0, return; end
            
            % 构造测量向量 y_meas 与 R 协方差
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
            R_cov = diag(R_list);
            
            % --- 2. 在当前先验预测状态下，单次计算残差与雅可比 (无需迭代) [4] ---
            h_val = zeros(M, 1);
            H_jac = zeros(M, 12*I); % 12I 维的局部测量流形误差雅可比
            
            row_ptr = 1;
            
            % A. 加速度计雅可比
            for r = 1:M_acc
                i = active_acc(r, 1);
                R_i = obj.states(i).R;
                a_i = obj.states(i).a;
                
                h_val(row_ptr:row_ptr+2) = R_i' * (a_i - obj.g_vec);
                H_jac(row_ptr:row_ptr+2, (i-1)*12 + (4:6)) = R_i';
                H_jac(row_ptr:row_ptr+2, (i-1)*12 + (7:9)) = skew(R_i' * (a_i - obj.g_vec));
                row_ptr = row_ptr + 3;
            end
            
            % B. 陀螺仪雅可比
            for r = 1:M_gyro
                i = active_gyro(r, 1);
                omega_i = obj.states(i).omega;
                
                h_val(row_ptr:row_ptr+2) = omega_i;
                H_jac(row_ptr:row_ptr+2, (i-1)*12 + (10:12)) = eye(3);
                row_ptr = row_ptr + 3;
            end
            
            % C. UWB 基站测距雅可比
            for r = 1:M_anc
                i = active_anc(r, 1); k = active_anc(r, 2);
                p_i = obj.states(i).p;
                c_k = anchors(k, :)';
                
                dist = norm(p_i - c_k);
                if dist < 1e-6, dist = 1e-6; end
                h_val(row_ptr) = dist;
                
                u_dir = (p_i - c_k)' / dist;
                H_jac(row_ptr, (i-1)*12 + (1:3)) = u_dir;
                row_ptr = row_ptr + 1;
            end
            
            % D. UWB 车间测距雅可比
            for r = 1:M_rel
                i = active_rel(r, 1); j = active_rel(r, 2);
                p_i = obj.states(i).p;
                p_j = obj.states(j).p;
                
                dist = norm(p_i - p_j);
                if dist < 1e-6, dist = 1e-6; end
                h_val(row_ptr) = dist;
                
                u_dir = (p_i - p_j)' / dist;
                H_jac(row_ptr, (i-1)*12 + (1:3)) = u_dir;
                H_jac(row_ptr, (j-1)*12 + (1:3)) = -u_dir;
                row_ptr = row_ptr + 1;
            end
            
            % 测量残差 (观测创新)
            y_err = y_meas - h_val;
            
            % --- 3. 维度转换投影 (将 12I 维雅可比映射至全维 15I 切空间) [7] ---
            % 构造状态选择矩阵 \pi (12I x 15I)
            pi_mat = zeros(12*I, 15*I);
            for i = 1:I
                pi_i = zeros(12, 15);
                pi_i(1:3, 1:3)     = eye(3);  % \delta p
                pi_i(4:6, 7:9)     = eye(3);  % \delta a
                pi_i(7:9, 10:12)   = eye(3);  % \delta \phi
                pi_i(10:12, 13:15) = eye(3);  % \delta \omega
                
                row_idx = (i-1)*12 + (1:12);
                col_idx = (i-1)*15 + (1:15);
                pi_mat(row_idx, col_idx) = pi_i;
            end
            
            % 全维 15I 状态量对应的测量雅可比 [7]
            H_full = H_jac * pi_mat; 
            
            % --- 4. 标准卡尔曼形式更新 ---
            S = H_full * obj.P * H_full' + R_cov;
            K = (obj.P * H_full') / S;
            
            % 求解 15I 维状态更新偏差 dx
            dx = K * y_err;
            
            % 更新系统状态协方差 P (带对称性保护)
            obj.P = (eye(15*I) - K * H_full) * obj.P;
            obj.P = 0.5 * (obj.P + obj.P');
            
            % --- 5. 标称状态流形回馈修正 [7] ---
            for i = 1:I
                dx_i = dx((i-1)*15 + (1:15));
                obj.states(i) = manifold_add(obj.states(i), dx_i);
            end
        end
    end
end