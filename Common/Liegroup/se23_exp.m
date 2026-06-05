function X = se23_exp(xi)
    % SE23_EXP 计算 SE2(3) 的指数映射
    % 输入: xi - 9x1 切向量 [phi; dv; dp]，其中各分量为 3x1 向量
    % 输出: X - 5x5 的 SE2(3) 齐次变换矩阵 (对应文档 Eq 5)
    
    phi = xi(1:3);
    dv  = xi(4:6);
    dp  = xi(7:9);
    
    theta = norm(phi);
    phi_hat = skew(phi);
    
    % Rodrigues 公式计算 SO(3) 旋转矩阵 R = exp(phi_hat)
    if theta < 1e-6
        R = eye(3) + phi_hat;
    else
        R = eye(3) + (sin(theta) / theta) * phi_hat + ...
            ((1 - cos(theta)) / (theta^2)) * (phi_hat^2);
    end
    
    % 计算 SO(3) 左雅可比
    Jl = so3_left_jacobian(phi);
    
    % 构造 5x5 矩阵
    X = eye(5);
    X(1:3, 1:3) = R;
    X(1:3, 4)   = Jl * dv;
    X(1:3, 5)   = Jl * dp;
end