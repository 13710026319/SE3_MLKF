function R = so3_exp(phi)
    % SO3_EXP SO(3)群的指数映射 (Rodrigues公式)
    % 输入: phi - 3x1 旋转向量 (轴角)
    % 输出: R - 3x3 SO(3) 旋转矩阵 (对应文档 Eq 12)
    
    theta = norm(phi);
    phi_skew = skew(phi);
    
    if theta < 1e-6
        % 极小值下的泰勒展开近似
        R = eye(3) + phi_skew;
    else
        R = eye(3) + (sin(theta) / theta) * phi_skew + ...
            ((1 - cos(theta)) / (theta^2)) * (phi_skew^2);
    end
end