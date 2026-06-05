classdef ESEKF < manifold_predict
    methods
        function obj = ESEKF(init_states_struct, init_P, IMU_noise_params)
            % 调用父类构造函数初始化状态
            obj@manifold_predict(init_states_struct, init_P, IMU_noise_params);
        end
        
        function update(obj, anchors, anc_meas, rel_meas, sigma_anc, sigma_rel)
            % ESEKF 观测更新步骤
            % 输入:
            %   anchors:   [Anchor_num x 3] 3D 基站位置矩阵
            %   anc_meas:  [N x Anchor_num] 当前时刻各车对基站的测距矩阵 (NaN 表示无观测)
            %   rel_meas:  [N x N] 当前时刻车间相对测距矩阵 (NaN 表示无观测)
            %   sigma_anc: 基站测距噪声标准差
            %   sigma_rel: 相对测距噪声标准差
            
            N = obj.Vehicle_num;
            
            H_list = {};
            y_list = {};
            R_list = [];
            
            % --- A. 线性化基站测距观测 ---
            for n = 1:N
                p_n = obj.states(n).X(1:3, 5); % 标称位置 [1]
                R_n = obj.states(n).X(1:3, 1:3);
                
                for a_idx = 1:size(anchors, 1)
                    d_meas = anc_meas(n, a_idx);
                    if isnan(d_meas), continue; end % 无有效观测，跳过
                    
                    % 预测测距 [64]
                    d_pred = norm(p_n - anchors(a_idx, :)');
                    
                    % 观测残差
                    dy = d_meas - d_pred;
                    
                    % 计算雅可比：测距关于位置的导数 dp_dir
                    dp_dir = (p_n - anchors(a_idx, :)')' / d_pred; % [1 x 3] 航向矢量
                    
                    % 构建联合误差状态 (15N 维) 的 H 矩阵行向量 (根据文档 Eq 83 的位置旋转提取)
                    H_row = zeros(1, 15*N);
                    % 位置误差状态排在 15D 单车误差状态的第 7:9 维 [89]
                    H_row((n-1)*15 + (7:9)) = dp_dir * R_n;
                    
                    H_list{end+1} = H_row;
                    y_list{end+1} = dy;
                    R_list(end+1) = sigma_anc^2;
                end
            end
            
            % --- B. 线性化车间相对测距观测 (避免双向重复计数) ---
            for n = 1:N
                p_n = obj.states(n).X(1:3, 5);
                R_n = obj.states(n).X(1:3, 1:3);
                
                for m = (n+1):N
                    d_meas = rel_meas(n, m);
                    if isnan(d_meas), continue; end
                    
                    p_m = obj.states(m).X(1:3, 5);
                    R_m = obj.states(m).X(1:3, 1:3);
                    
                    % 预测相对距离 [65]
                    d_pred = norm(p_n - p_m);
                    dy = d_meas - d_pred;
                    
                    % 测距方向向量
                    dp_dir = (p_n - p_m)' / d_pred;
                    
                    H_row = zeros(1, 15*N);
                    % 填充车 n 的位置导数分量 (正方向) [73]
                    H_row((n-1)*15 + (7:9)) = dp_dir * R_n;
                    % 填充车 m 的位置导数分量 (反方向) [73]
                    H_row((m-1)*15 + (7:9)) = -dp_dir * R_m;
                    
                    H_list{end+1} = H_row;
                    y_list{end+1} = dy;
                    R_list(end+1) = sigma_rel^2;
                end
            end
            
            % 如果当前时刻没有任何传感器观测，直接返回
            if isempty(H_list), return; end
            
            % 拼接为标准的卡尔曼矩阵形式
            H = cell2mat(H_list');
            Y = cell2mat(y_list');
            R_cov = diag(R_list);
            
            % 3. 卡尔曼增益与后验计算
            S = H * obj.P * H' + R_cov;
            K = obj.P * H' / S;
            
            % 计算后验误差状态量
            dx = K * Y;
            
            % 更新协方差矩阵 P (采用对称正定化保护)
            I_joint = eye(15*N);
            obj.P = (I_joint - K * H) * obj.P;
            obj.P = 0.5 * (obj.P + obj.P');
            
            % 4. 复合流形状态修正并复位误差状态 (文档 4.3 节)
            for n = 1:N
                dx_n = dx((n-1)*15 + (1:15));
                
                % 调用此前编写的 manifold_add 流形加法函数 [11, 13]
                obj.states(n) = manifold_add(obj.states(n), dx_n);
            end
            
            % 更新完成后，后验误差状态已被重置（流形更新的反馈性质），直接进入下一个传播步骤 [11]。
        end
    end
end