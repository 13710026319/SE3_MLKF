function AdX = se23_adjoint(X)
    % SE23_ADJOINT 计算 SE2(3) 的伴随矩阵
    % 输入: X - 5x5 的 SE2(3) 齐次变换矩阵
    % 输出: AdX - 9x9 的伴随变换矩阵 (对应文档 Eq 8)
    
    R = X(1:3, 1:3);
    v = X(1:3, 4);
    p = X(1:3, 5);
    
    v_hat = skew(v);
    p_hat = skew(p);
    
    AdX = zeros(9, 9);
    AdX(1:3, 1:3) = R;
    AdX(4:6, 1:3) = v_hat * R;
    AdX(4:6, 4:6) = R;
    AdX(7:9, 1:3) = p_hat * R;
    AdX(7:9, 7:9) = R;
end