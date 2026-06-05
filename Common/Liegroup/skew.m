function phi_hat = skew(phi)
    % SKEW 计算三维向量的反对称矩阵
    % 输入: phi - 3x1 向量
    % 输出: phi_hat - 3x3 反对称矩阵
    
    phi_hat = [      0, -phi(3),  phi(2);
                phi(3),       0, -phi(1);
               -phi(2),  phi(1),       0];
end