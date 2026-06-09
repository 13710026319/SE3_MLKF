function Jr_inv = so3_inv_right_jacobian(phi)
    % SO3_INV_RIGHT_JACOBIAN 计算 SO(3) 的逆右雅可比矩阵
    % 输入: phi - 3x1 旋转向量
    % 输出: Jr_inv - 3x3 逆右雅可比矩阵 (对应文档 Eq 54)
    
    theta = norm(phi);
    phi_skew = skew(phi);
    
    if theta < 1e-6
        % 极小值下的极限近似值 (cot(x) 展开后抵消所得)
        coeff = 1/12;
    else
        coeff = 1 / (theta^2) - cot(theta/2) / (2 * theta);
    end
    
    Jr_inv = eye(3) + 0.5 * phi_skew + coeff * (phi_skew^2);
end