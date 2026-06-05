function xi = se23_log(X)
    % SE23_LOG 计算 SE2(3) 的对数映射
    % 输入: X - 5x5 的 SE2(3) 齐次变换矩阵
    % 输出: xi - 9x1 切向量 [phi; dv; dp] (对应文档 Eq 7)
    
    R = X(1:3, 1:3);
    v = X(1:3, 4);
    p = X(1:3, 5);
    
    % 计算旋转向量 phi
    tr = trace(R);
    val = (tr - 1) / 2;
    val = max(-1, min(1, val)); % 数值防越界安全截断
    theta = acos(val);
    
    R_diff = R - R';
    if theta < 1e-6
        phi = unskew(R_diff) / 2;
    else
        phi = (theta / (2 * sin(theta))) * unskew(R_diff);
    end
    
    % 计算 SO(3) 逆左雅可比 Jl_inv
    phi_hat = skew(phi);
    if theta < 1e-6
        Jl_inv = eye(3) - 0.5 * phi_hat + (1/12) * (phi_hat^2);
    else
        Jl_inv = eye(3) - 0.5 * phi_hat + ...
            ((1 / (theta^2)) - (1 + cos(theta)) / (2 * theta * sin(theta))) * (phi_hat^2);
    end
    
    % 还原平移与速度摄动量
    dv = Jl_inv * v;
    dp = Jl_inv * p;
    
    xi = [phi; dv; dp];
end