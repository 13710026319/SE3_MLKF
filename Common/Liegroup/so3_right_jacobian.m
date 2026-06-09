function Jr = so3_right_jacobian(phi)
    % SO3_RIGHT_JACOBIAN 计算 SO(3) 的标称右雅可比矩阵
    % 输入: phi - 3x1 旋转向量
    % 输出: Jr - 3x3 右雅可比矩阵 (对应文档 Eq 28)
    
    theta = norm(phi);
    phi_skew = skew(phi);
    
    if theta < 1e-6
        % 极小值下的极限近似 (罗德里格斯展开前两项)
        Jr = eye(3) - 0.5 * phi_skew + (1/6) * (phi_skew^2);
    else
        Jr = eye(3) - ((1 - cos(theta)) / (theta^2)) * phi_skew + ...
             ((theta - sin(theta)) / (theta^3)) * (phi_skew^2);
    end
end