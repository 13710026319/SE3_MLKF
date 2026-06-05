function Jl = so3_left_jacobian(phi)
    % SO3_LEFT_JACOBIAN 计算 SO(3) 的左雅可比矩阵
    % 输入: phi - 3x1 旋转向量
    % 输出: Jl - 3x3 左雅可比矩阵 (对应文档 Eq 6)
    
    theta = norm(phi);
    phi_hat = skew(phi);
    
    if theta < 1e-6
        % 极小值下的极限近似
        Jl = eye(3) + 0.5 * phi_hat;
    else
        Jl = eye(3) + ((1 - cos(theta)) / (theta^2)) * phi_hat + ...
             ((theta - sin(theta)) / (theta^3)) * (phi_hat^2);
    end
end