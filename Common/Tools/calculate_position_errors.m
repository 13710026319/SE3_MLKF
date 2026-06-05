function [errors, rmse] = calculate_position_errors(pos_est, pos_true)
    % CALCULATE_POSITION_ERRORS 计算估计位置和真实位置之间的三轴误差、欧氏距离误差及RMSE
    % 输入:
    %   pos_est:  Cell 数组 {Vehicle_num x 1}，每个单元为 [N_steps x 3] 的估计位置矩阵
    %   pos_true: Cell 数组 {Vehicle_num x 1}，每个单元为 [N_steps x 3] 的真实位置矩阵
    % 输出:
    %   errors:   结构体数组 (1 x Vehicle_num)，每个元素包含:
    %       .axis_err: [N_steps x 3] 的三轴位置误差 (X, Y, Z)
    %       .euc_err:  [N_steps x 1] 的欧氏距离误差
    %   rmse:     结构体数组 (1 x Vehicle_num)，每个元素包含:
    %       .axis_rmse: [1 x 3] 三轴位置的 RMSE
    %       .euc_rmse:  [1 x 1] 欧氏距离的 RMSE

    Vehicle_num = length(pos_est);
    
    % 初始化输出结构体
    errors = struct('axis_err', {}, 'euc_err', {});
    rmse = struct('axis_rmse', {}, 'euc_rmse', {});
    
    for n = 1:Vehicle_num
        % 1. 计算三轴位置误差 (Estimated - True)
        axis_err = pos_est{n} - pos_true{n};
        
        % 2. 计算三轴平方距离之和并开方，获得欧氏定位误差
        euc_err = sqrt(sum(axis_err.^2, 2));
        
        % 写入误差曲线序列
        errors(n).axis_err = axis_err;
        errors(n).euc_err = euc_err;
        
        % 3. 计算三轴方向的 RMSE [2]
        axis_rmse = sqrt(mean(axis_err.^2, 1));
        
        % 4. 计算欧氏距离的 RMSE [2]
        euc_rmse = sqrt(mean(euc_err.^2));
        
        % 写入结果
        rmse(n).axis_rmse = axis_rmse;
        rmse(n).euc_rmse = euc_rmse;
    end
end