classdef manifold_predict < handle
    properties
        Vehicle_num        % 车辆数量
        states             % 结构体数组：states(n).X (5x5 李群), states(n).b (6x1 零偏)
        P                  % 联合误差状态协方差矩阵 [15N x 15N]
        Q_imu_params       % IMU 噪声参数结构体
    end
    
    methods
        function obj = manifold_predict(init_states_struct, init_P, IMU_noise_params)
            % 构造函数
            obj.Vehicle_num = length(init_states_struct);
            obj.states = init_states_struct;
            obj.P = init_P;
            obj.Q_imu_params = IMU_noise_params;
        end
        
        function propagate(obj, imu_acc, imu_gyro, dt)
            % 状态与协方差传播 (对应文档 3.4 节)
            
            N = obj.Vehicle_num;
            g_vec = [0; 0; -9.81]; % 3D重力矢量 (3x1 向量)
            
            % 预分配联合状态转移和过程噪声矩阵
            Phi_joint = zeros(15*N, 15*N);
            Q_joint = zeros(15*N, 15*N);
            
            for n = 1:N
                % 1. 获取输入与标称零偏
                acc_m = imu_acc(n, :)';
                gyro_m = imu_gyro(n, :)';
                
                ba = obj.states(n).b(1:3);
                bg = obj.states(n).b(4:6);
                
                % 2. 标称零偏修正 (文档 Eq 29-30)
                omega_hat = gyro_m - bg;
                a_hat = acc_m - ba;
                
                % 3. 标称状态时间传播 (文档 Eq 53-57)
                R_prev = obj.states(n).X(1:3, 1:3);
                v_prev = obj.states(n).X(1:3, 4);
                p_prev = obj.states(n).X(1:3, 5);
                
                % 【修正部分】: 旋转传播，直接计算 SO(3) 上的罗德里格斯指数映射 (3x3 矩阵) [3]
                phi_step = omega_hat * dt;
                theta_omega = norm(phi_step);
                omega_skew = skew(phi_step);
                
                if theta_omega < 1e-6
                    exp_R = eye(3) + omega_skew;
                else
                    exp_R = eye(3) + (sin(theta_omega) / theta_omega) * omega_skew + ...
                            ((1 - cos(theta_omega)) / (theta_omega^2)) * (omega_skew^2);
                end
                
                % 正确执行 3x3 矩阵相乘 [7]
                R_next = R_prev * exp_R;
                
                % 速度传播 [54]
                v_next = v_prev + (R_prev * a_hat + g_vec) * dt;
                
                % 位置传播 [55]
                p_next = p_prev + v_prev * dt + 0.5 * (R_prev * a_hat + g_vec) * dt^2;
                
                % 更新标称状态
                obj.states(n).X(1:3, 1:3) = R_next;
                obj.states(n).X(1:3, 4)   = v_next;
                obj.states(n).X(1:3, 5)   = p_next;
                
                % 4. 计算该车的离散状态转移矩阵 Phi_n (文档 Eq 59)
                omega_skew_meas = skew(omega_hat);
                a_skew = skew(a_hat);
                I3 = eye(3);
                
                Phi_n = eye(15);
                % 第一行
                Phi_n(1:3, 1:3)   = I3 - omega_skew_meas * dt;
                Phi_n(1:3, 13:15) = -I3 * dt;
                % 第二行
                Phi_n(4:6, 1:3)   = -a_skew * dt;
                Phi_n(4:6, 4:6)   = I3 - omega_skew_meas * dt;
                Phi_n(4:6, 10:12) = -I3 * dt;
                % 第三行
                Phi_n(7:9, 4:6)   = I3 * dt;
                Phi_n(7:9, 7:9)   = I3 - omega_skew_meas * dt;
                
                % 5. 计算该车的离散过程噪声协方差 Q_dn (文档 Eq 52, 62)
                sig_na = obj.Q_imu_params.sigma_na;
                sig_nw = obj.Q_imu_params.sigma_nw;
                sig_ba = obj.Q_imu_params.sigma_ba;
                sig_bw = obj.Q_imu_params.sigma_bw;
                
                Q_dn = zeros(15, 15);
                Q_dn(1:3, 1:3)     = (sig_nw^2 * dt) * I3;   % 角度噪声
                Q_dn(4:6, 4:6)     = (sig_na^2 * dt) * I3;   % 速度噪声
                Q_dn(10:12, 10:12) = (sig_ba^2 * dt) * I3;   % 加速度零偏随机游走
                Q_dn(13:15, 13:15) = (sig_bw^2 * dt) * I3;   % 陀螺仪零偏随机游走
                
                % 填入对角块
                idx = (n-1)*15 + (1:15);
                Phi_joint(idx, idx) = Phi_n;
                Q_joint(idx, idx) = Q_dn;
            end
            
            % 6. 联合协方差传播 (文档 Eq 63)
            obj.P = Phi_joint * obj.P * Phi_joint' + Q_joint;
            
            % 保持对称性
            obj.P = 0.5 * (obj.P + obj.P');
        end
    end
end
