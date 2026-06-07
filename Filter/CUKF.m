classdef CUKF < manifold_predict
    properties
        alpha   % UKF 比例因子 (通常设为 1e-3 ~ 1)
        beta    % 状态分布参数 (对于高斯分布，常用 2)
        kappa   % 辅助比例因子 (常用 0)
    end
    
    methods
        function obj = CUKF(init_states_struct, init_P, IMU_noise_params)
            % 调用父类 manifold_predict 构造函数
            obj@manifold_predict(init_states_struct, init_P, IMU_noise_params);
            
            % 设置默认的无迹变换参数
            obj.alpha = 1e-3; 
            obj.beta = 2;
            obj.kappa = 0;
        end
        
        function update(obj, anchors, anc_meas, rel_meas, sigma_anc, sigma_rel)
            % CUKF 集中式无迹卡尔曼更新步 (在 15N 维切空间撒点)
            
            N = obj.Vehicle_num;
            L = 15 * N; % 误差状态联合维度 (如果是4辆车，则为 60 维)
            
            % 1. 整理并统计活跃测距链路
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
            
            if isempty(active_anchors) && isempty(active_rel)
                return;
            end
            
            % 组装测量向量
            y_meas = [active_anchors(:, 3); active_rel(:, 3)];
            M_anc = size(active_anchors, 1);
            M_rel = size(active_rel, 1);
            M = M_anc + M_rel;
            
            R_list = [ones(M_anc, 1) * sigma_anc^2; ones(M_rel, 1) * sigma_rel^2];
            R_cov = diag(R_list); % 测量噪声协方差矩阵
            
            % 2. 计算无迹变换参数和权重
            lambda = obj.alpha^2 * (L + obj.kappa) - L;
            
            Wm = zeros(2*L+1, 1);
            Wc = zeros(2*L+1, 1);
            
            Wm(1) = lambda / (L + lambda);
            Wc(1) = lambda / (L + lambda) + (1 - obj.alpha^2 + obj.beta);
            
            common_weight = 1 / (2 * (L + lambda));
            Wm(2:end) = common_weight;
            Wc(2:end) = common_weight;
            
            % 3. 生成切空间(误差状态) Sigma 点，并利用 sqrtm/SVD 抵御奇异不确定度
            try
                % Cholesky 分解
                chol_P = chol((L + lambda) * obj.P);
                S = chol_P'; % 获得下三角，每一列代表一个扰动方向
            catch
                % 若协方差由于不可观测方向导致半正定，采用数值更强的 SVD 求解平方根
                [U, S_val, ~] = svd((L + lambda) * obj.P);
                S = U * sqrt(max(0, S_val));
            end
            
            % 在 15N 维切平面生成 2L+1 个误差状态 Sigma 向量
            chi_points = zeros(L, 2*L+1);
            chi_points(:, 1) = 0;
            chi_points(:, 2:L+1) = S;
            chi_points(:, L+2:2*L+1) = -S;
            
            % 4. 将切空间扰动映射到李群流形，生成流形状态 Sigma 集合 [4, 13]
            x_sigma = cell(2*L+1, 1);
            for i = 1:(2*L+1)
                dx_joint = chi_points(:, i);
                
                % 对联合状态的每辆车，单独进行流形广义加法 (x_nominal ⊕ dx) [13]
                x_sigma{i} = struct('X', {}, 'b', {});
                for n = 1:N
                    dx_n = dx_joint((n-1)*15 + (1:15));
                    x_sigma{i}(n) = manifold_add(obj.states(n), dx_n);
                end
            end
            
            % 5. 非线性测量变换：将流形 Sigma 状态通过观测函数 h(x) [8]
            gamma = zeros(M, 2*L+1);
            for i = 1:(2*L+1)
                curr_state = x_sigma{i};
                idx_m = 1;
                
                % 基站测距
                for r = 1:M_anc
                    n = active_anchors(r, 1);
                    anc_idx = active_anchors(r, 2);
                    p_n = curr_state(n).X(1:3, 5);
                    gamma(idx_m, i) = norm(p_n - anchors(anc_idx, :)');
                    idx_m = idx_m + 1;
                end
                
                % 车间测距
                for r = 1:M_rel
                    n = active_rel(r, 1);
                    m = active_rel(r, 2);
                    p_n = curr_state(n).X(1:3, 5);
                    p_m = curr_state(m).X(1:3, 5);
                    gamma(idx_m, i) = norm(p_n - p_m);
                    idx_m = idx_m + 1;
                end
            end
            
            % 6. 计算均值与协方差
            % 预测测量均值 y_pred
            y_pred = zeros(M, 1);
            for i = 1:(2*L+1)
                y_pred = y_pred + Wm(i) * gamma(:, i);
            end
            
            % 预测测量协方差 Pyy 与 互协方差 Pxy
            Pyy = zeros(M, M);
            Pxy = zeros(L, M);
            for i = 1:(2*L+1)
                dy = gamma(:, i) - y_pred;
                Pyy = Pyy + Wc(i) * (dy * dy');
                
                % 切平面误差状态扰动项即为交叉项 (chi_points_i - 0)
                dx_tangent = chi_points(:, i);
                Pxy = Pxy + Wc(i) * (dx_tangent * dy');
            end
            
            % 叠加测量噪声
            Pyy = Pyy + R_cov;
            
            % 7. 卡尔曼增益与后验状态更新
            K = Pxy / Pyy;
            
            % 计算切空间修正增量 dx
            dx = K * (y_meas - y_pred);
            
            % 更新联合协方差 P (Joseph 等价稳定性更新形式，防止非正定)
            obj.P = obj.P - K * Pyy * K';
            obj.P = 0.5 * (obj.P + obj.P');
            
            % 8. 标称状态流形纠正 [11, 13]
            for n = 1:N
                dx_n = dx((n-1)*15 + (1:15));
                obj.states(n) = manifold_add(obj.states(n), dx_n);
            end
        end
    end
end