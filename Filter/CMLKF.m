classdef CMLKF < manifold_predict
    methods
        function obj = CMLKF(init_states_struct, init_P, IMU_noise_params)
            % 调用父类 manifold_predict 构造函数
            obj@manifold_predict(init_states_struct, init_P, IMU_noise_params);
        end
        
        function update(obj, anchors, anc_meas, rel_meas, sigma_anc, sigma_rel)
            
            N = obj.Vehicle_num;
            
            % 1. 整理当前时刻所有车辆活跃的测距链路
            active_anchors = []; % [vehicle_idx, anchor_idx, distance_val]
            for n = 1:N
                for i = 1:size(anchors, 1)
                    if ~isnan(anc_meas(n, i))
                        active_anchors = [active_anchors; n, i, anc_meas(n, i)];
                    end
                end
            end

            active_rel = []; % [veh_n, veh_m, distance_val]
            for n = 1:N
                for m = (n+1):N
                    if ~isnan(rel_meas(n, m))
                        active_rel = [active_rel; n, m, rel_meas(n, m)];
                    end
                end
            end

            % 如果当前时刻无任何有效观测，直接返回
            if isempty(active_anchors) && isempty(active_rel)
                return;
            end
            

            % 2. 构造联合测量向量 y_meas 与 协方差逆矩阵 R_inv
            y_meas = [active_anchors(:, 3); active_rel(:, 3)];
            M_anc = size(active_anchors, 1);
            M_rel = size(active_rel, 1);
            M = M_anc + M_rel;
            
            R_list = [ones(M_anc, 1) * sigma_anc^2; ones(M_rel, 1) * sigma_rel^2];
            R_inv = diag(1 ./ R_list); % 测量逆协方差矩阵 (对角线形式加速运算)
            
            % 3. 联合位置优化 (Levenberg-Marquardt 优化形式)
            % 标称位置预测值作为初值 p^(0) (Eq 69)
            p_pred = zeros(3*N, 1);
            for n = 1:N
                p_pred((n-1)*3 + (1:3)) = obj.states(n).X(1:3, 5);
            end
            p_opt = p_pred; 
            
            max_iter = 10;
            tol = 1e-4;
            
            for iter = 1:max_iter
                h_val = zeros(M, 1);
                H_jac = zeros(M, 3*N);
                
                % A. 计算基站观测模型及其对三轴位置的导数
                for r_idx = 1:M_anc
                    n = active_anchors(r_idx, 1);
                    i = active_anchors(r_idx, 2);
                    
                    p_n = p_opt((n-1)*3 + (1:3));
                    a_i = anchors(i, :)';
                    
                    dist = norm(p_n - a_i);
                    if dist < 1e-6, dist = 1e-6; end
                    
                    h_val(r_idx) = dist;
                    dp_dir = (p_n - a_i)' / dist;
                    H_jac(r_idx, (n-1)*3 + (1:3)) = dp_dir;
                end
                
                % B. 计算车间观测模型及其对三轴位置的导数
                for r_idx = 1:M_rel
                    n = active_rel(r_idx, 1);
                    m = active_rel(r_idx, 2);
                    
                    p_n = p_opt((n-1)*3 + (1:3));
                    p_m = p_opt((m-1)*3 + (1:3));
                    
                    dist = norm(p_n - p_m);
                    if dist < 1e-6, dist = 1e-6; end
                    
                    h_val(M_anc + r_idx) = dist;
                    dp_dir = (p_n - p_m)' / dist;
                    H_jac(M_anc + r_idx, (n-1)*3 + (1:3)) = dp_dir;
                    H_jac(M_anc + r_idx, (m-1)*3 + (1:3)) = -dp_dir;
                end
                
                % 测量残差 (Eq 75)
                residual = y_meas - h_val;
                
                % 计算位置增量 (Eq 78)
                H_T_R_inv = H_jac' * R_inv;
                Hessian_pos = H_T_R_inv * H_jac;
                
                % 【核心改进 1】：引入 Levenberg-Marquardt 阻尼因子 (1e-4)
                % 保证即使在高度维度、2基站完全不可观测时，Hessian 也严格正定，消除所有奇异和 NaN 警告
                Hessian_pos = Hessian_pos + 1e-4 * eye(3*N);
                
                Delta_p = Hessian_pos \ (H_T_R_inv * residual);
                p_opt = p_opt + Delta_p;
                
                if norm(Delta_p) < tol
                    break;
                end
            end
            
            % GN 优化收敛，记录 ML 估计最优结果
            p_ML = p_opt;
            
            % 【修正一】：在最终收敛点 p_ML 处，重新严格计算一次雅可比 H_jac [10]
            H_jac_ML = zeros(M, 3*N);
            for r_idx = 1:M_anc
                n = active_anchors(r_idx, 1);
                i = active_anchors(r_idx, 2);
                p_n = p_ML((n-1)*3 + (1:3));
                a_i = anchors(i, :)';
                dist = norm(p_n - a_i);
                if dist < 1e-6, dist = 1e-6; end
                H_jac_ML(r_idx, (n-1)*3 + (1:3)) = (p_n - a_i)' / dist;
            end
            for r_idx = 1:M_rel
                n = active_rel(r_idx, 1);
                m = active_rel(r_idx, 2);
                p_n = p_ML((n-1)*3 + (1:3));
                p_m = p_ML((m-1)*3 + (1:3));
                dist = norm(p_n - p_m);
                if dist < 1e-6, dist = 1e-6; end
                dp_dir = (p_n - p_m)' / dist;
                H_jac_ML(M_anc + r_idx, (n-1)*3 + (1:3)) = dp_dir;
                H_jac_ML(M_anc + r_idx, (m-1)*3 + (1:3)) = -dp_dir;
            end
            
            % 计算严格对应 p_ML 的收敛 Hessian 矩阵 (Eq 80) [10]
            H_T_R_inv = H_jac_ML' * R_inv;
            Xi_k = H_T_R_inv * H_jac_ML;
            
            % 4. 计算位置误差投影矩阵 D_k (Eq 82-83)
            D_k = zeros(3*N, 15*N);
            for n = 1:N
                R_n = obj.states(n).X(1:3, 1:3);
                row_idx = (n-1)*3 + (1:3);
                col_idx = (n-1)*15 + (7:9);
                D_k(row_idx, col_idx) = R_n;
            end
            
            % 【核心改进 2】：对虚拟测量噪声 Xi_k 施加正定性保护 (加上 1e-4 阻尼)
            % 在高度、基站不可观测维度上，相当于配置了一个高噪声虚拟观测（方差为10^4），从而保持数值绝对稳定 [10]
            R_eff = (Xi_k + 1e-4 * eye(3*N)) \ eye(3*N); 
            
            % 等效卡尔曼滤波协方差形式更新 [10]
            S_eff = D_k * obj.P * D_k' + R_eff;
            K_eff = (obj.P * D_k') / S_eff;
            
            % 虚拟观测残差 (Eq 86)
            delta_y_k = p_ML - p_pred;
            
            % 后验状态偏差 dx 更新 [10]
            dx = K_eff * delta_y_k;
            
            % 后验联合协方差 P 更新 [10]
            obj.P = (eye(15*N) - K_eff * D_k) * obj.P;
            obj.P = 0.5 * (obj.P + obj.P');
            
            % 6. 流形标称状态纠正并重置误差 [11, 13]
            for n = 1:N
                dx_n = dx((n-1)*15 + (1:15));
                obj.states(n) = manifold_add(obj.states(n), dx_n);
            end
        end
    end
end