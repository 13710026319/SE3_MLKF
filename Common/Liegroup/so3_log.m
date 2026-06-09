function phi = so3_log(R)
    % SO3_LOG SO(3)群的对数映射
    % 输入: R - 3x3 SO(3) 旋转矩阵
    % 输出: phi - 3x1 旋转向量
    
    tr = trace(R);
    val = (tr - 1) / 2;
    val = max(-1, min(1, val)); % 数值防越界截断
    theta = acos(val);
    
    R_diff = R - R';
    if theta < 1e-6
        phi = unskew(R_diff) / 2;
    else
        phi = (theta / (2 * sin(theta))) * unskew(R_diff);
    end
end