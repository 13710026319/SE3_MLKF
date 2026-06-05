function phi = unskew(phi_hat)
    % UNSKEW 从 3x3 反对称矩阵中提取三维向量
    % 输入: phi_hat - 3x3 反对称矩阵
    % 输出: phi - 3x1 向量
    
    phi = [phi_hat(3, 2); 
           phi_hat(1, 3); 
           phi_hat(2, 1)];
end