function dx = manifold_sub(x1, x2)
    % MANIFOLD_SUB 复合流形 M = SE2(3) x R^6 上的广义减法 (Minus)
    % 计算表达式为: dx = x1 ⊖ x2 (对应文档 Eq 14)
    % 输入:
    %   x1, x2 - 两个状态结构体 (包含 .X 和 .b)
    % 输出:
    %   dx     - 15x1 的流形误差向量 [dchi; db]
    
    % 提取并解析 X2
    R2 = x2.X(1:3, 1:3);
    v2 = x2.X(1:3, 4);
    p2 = x2.X(1:3, 5);
    
    % 解析法求 X2_inv (解析求逆提升运算鲁棒性)
    X2_inv = eye(5);
    X2_inv(1:3, 1:3) = R2';
    X2_inv(1:3, 4)   = -R2' * v2;
    X2_inv(1:3, 5)   = -R2' * p2;
    
    % 计算李群差值 X_diff = X2_inv * X1
    X_diff = X2_inv * x1.X;
    
    % 投影到切空间
    dchi = se23_log(X_diff);
    
    % 偏置线性相减
    db = x1.b - x2.b;
    
    dx = [dchi; db];
end