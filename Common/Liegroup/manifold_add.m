function x_new = manifold_add(x_nominal, dx)
    % MANIFOLD_ADD 复合流形 M = SE2(3) x R^6 上的广义加法 (Plus)
    % 输入: 
    %   x_nominal - 标称状态结构体，包含：
    %       .X: 5x5 的 SE2(3) 矩阵
    %       .b: 6x1 的 IMU 零偏向量 [ba; bg]
    %   dx        - 15x1 的误差状态向量 [dchi; db] (dchi 为 9x1，db 为 6x1)
    % 输出:
    %   x_new     - 更新后的状态结构体 (对应文档 Eq 13)
    
    dchi = dx(1:9);
    db   = dx(10:15);
    
    % 更新李群部分：X_new = X * exp(dchi^)
    X_new = x_nominal.X * se23_exp(dchi);
    
    % 更新加性偏置部分
    b_new = x_nominal.b + db;
    
    x_new.X = X_new;
    x_new.b = b_new;
end