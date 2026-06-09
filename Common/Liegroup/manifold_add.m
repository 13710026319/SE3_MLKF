function x_new = manifold_add(x_nominal, dx)
    % MANIFOLD_ADD 复合状态流形 M 上的广义加法 (Plus) 运算
    % 输入: 
    %   x_nominal - 标称状态结构体，包含：
    %       .p: 3x1 位置向量 [1]
    %       .v: 3x1 速度向量 [1]
    %       .a: 3x1 加速度向量 [1]
    %       .R: 3x3 SO(3) 旋转矩阵 [1]
    %       .omega: 3x1 角速度向量 [1]
    %   dx        - 15x1 误差状态向量 [dp; dv; da; dphi; domega] [2]
    % 输出:
    %   x_new     - 更新后的复合标称状态结构体 (对应文档 Eq 13)
    
    dp     = dx(1:3);
    dv     = dx(4:6);
    da     = dx(7:9);
    dphi   = dx(10:12);
    domega = dx(13:15);
    
    % 分量执行流形直和更新
    x_new.p     = x_nominal.p + dp;
    x_new.v     = x_nominal.v + dv;
    x_new.a     = x_nominal.a + da;
    x_new.R     = x_nominal.R * so3_exp(dphi); % 旋转矩阵右乘指数映射进行Retraction [2]
    x_new.omega = x_nominal.omega + domega;
end