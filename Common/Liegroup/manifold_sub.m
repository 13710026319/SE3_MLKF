function dx = manifold_sub(x1, x2)
    % MANIFOLD_SUB 复合状态流形 M 上的广义减法 (Minus) 运算 (dx = x1 ⊟ x2)
    % 输入:
    %   x1, x2 - 两个复合状态结构体 (各自包含 .p, .v, .a, .R, .omega 分量)
    % 输出:
    %   dx     - 15x1 的切平面误差向量 [dp; dv; da; dphi; domega] (对应文档 Eq 14)
    
    dp     = x1.p - x2.p;
    dv     = x1.v - x2.v;
    da     = x1.a - x2.a;
    dphi   = so3_log(x2.R' * x1.R); % 旋转矩阵相减 (对应 log((R2)^T * R1)) [2]
    domega = x1.omega - x2.omega;
    
    dx = [dp; dv; da; dphi; domega];
end