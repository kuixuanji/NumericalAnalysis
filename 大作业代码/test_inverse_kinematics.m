function diagnostic_ik_test()
    clc; clear; close all;
    %% =====================================================================
    %% 1. 配置机械臂的 MDH 参数 (基于你的 8 行参数表)
    %% =====================================================================
    global MDH_params;
    MDH_params = [
        0,      pi,      -0.1564-0.1284;  % Base -> 轴1
        0,      pi/2,    -0.0054-0.0064;  % 轴1  -> 轴2
        0,     -pi/2,    -0.2104-0.2104;  % 轴2  -> 轴3
        0,      pi/2,    -0.0064-0.0064;  % 轴3  -> 轴4
        0,     -pi/2,    -0.2084-0.1059;  % 轴4  -> 轴5
        0,      pi/2,     0;              % 轴5  -> 轴6
        0,     -pi/2,    -0.0615-0.1059;  % 轴6  -> 轴7
        0,      pi,       0               % 轴7  -> Tool (第8行)
    ];
    %% =====================================================================
    %% 2. 设定测试的关节序列 (Ground Truth) 并生成目标 T 矩阵
    %% =====================================================================
    q_truth = [1, -0.3, 0.8, 1.2, 0.5, 0.6, 0]'; 
    T_target = Forward_Kinematics(q_truth);
    
    fprintf('==================================================\n');
    fprintf('           IK 独立分量误差曲线诊断程序            \n');
    fprintf('==================================================\n');
    
    %% =====================================================================
    %% 3. 运行带全数据记录的逆运动学迭代
    %% =====================================================================
    q_init = zeros(7, 1); % 初始猜测全 0
    [q_sol, success, err_history] = IK_Numerical_Diagnostic(T_target, q_init, MDH_params);
    %% =====================================================================
    %% 4. 绘制精细化的误差曲线图
    %% =====================================================================
    figure('Name', '逆运动学迭代误差分量诊断', 'NumberTitle', 'off', 'Position', [100, 100, 1100, 700]);
    iters = 1:size(err_history, 1);
    
    % 子图 1: 三维位置误差分量 (X, Y, Z) - 【单位已改为 mm】
    subplot(2, 2, 1);
    plot(iters, err_history(:, 1), '-r.', 'LineWidth', 1.5, 'DisplayName', 'e\_x'); hold on;
    plot(iters, err_history(:, 2), '-g.', 'LineWidth', 1.5, 'DisplayName', 'e\_y');
    plot(iters, err_history(:, 3), '-b.', 'LineWidth', 1.5, 'DisplayName', 'e\_z');
    grid on; xlabel('Iteration'); ylabel('Position Error (mm)');
    title('位置误差分量 (Position Error Components)');
    legend('Location', 'best');
    
    % 子图 2: 三维姿态误差分量 (wx, wy, wz)
    subplot(2, 2, 2);
    plot(iters, err_history(:, 4), '--r.', 'LineWidth', 1.5, 'DisplayName', 'e\_wx'); hold on;
    plot(iters, err_history(:, 5), '--g.', 'LineWidth', 1.5, 'DisplayName', 'e\_wy');
    plot(iters, err_history(:, 6), '--b.', 'LineWidth', 1.5, 'DisplayName', 'e\_wz');
    grid on; xlabel('Iteration'); ylabel('Orientation Error (rad)');
    title('姿态误差分量 (Orientation Error Components)');
    legend('Location', 'best');
    
    % 子图 3: 总误差模长收敛曲线 (Log 坐标) - 【位置部分以 mm 融合计算】
    subplot(2, 2, 3);
    total_norm = zeros(length(iters), 1);
    for k = 1:length(iters)
        total_norm(k) = norm(err_history(k, :));
    end
    semilogy(iters, total_norm, '-k^', 'LineWidth', 2, 'MarkerSize', 4);
    grid on; xlabel('Iteration'); ylabel('Total Error Norm (Mixed Units: mm & rad)');
    title('综合位姿误差模长收敛曲线');
    
    % 打印状态报告
    if success
        % 由于 err_history 内部前三列已经换算成 mm，这里直接取模即可
        fprintf('[🎉 成功] 算法收敛！最终位置残差: %.6f mm\n', norm(err_history(end, 1:3)));
    else
        fprintf('[❌ 失败] 算法未收敛。请观察右侧弹出的误差分量曲线图进行诊断。\n');
    end
end

%% =====================================================================
%% 诊断版逆运动学函数：实时记录 6 维误差数据 (位置已转换为 mm)
%% =====================================================================
function [q, tool, err_history] = IK_Numerical_Diagnostic(T_target, q_init, MDH_params)
    p_target = T_target(1:3, 4);       
    R_target = T_target(1:3, 1:3);     
    q_target_quat = rotMat2Quat(R_target); 
    
    q = q_init;
    max_iter = 50;
    tol = 1e-4;
    lambda_max = 0.02;
    tool = false;
    err_history = zeros(max_iter, 6); % 预分配空间存储 [ex, ey, ez, ewx, ewy, ewz]
    
    for iter = 1:max_iter
        T_curr = Forward_Kinematics(q);
        p_curr = T_curr(1:3, 4);
        R_curr = T_curr(1:3, 1:3);
        q_curr = rotMat2Quat(R_curr);
        
        % 1. 在算法求解中使用标准单位（米）以确保雅可比矩阵物理意义正确
        e_pos = p_target - p_curr;
        e_ori = Compute_QuatError(q_target_quat, q_curr);
        
        % 2. 【核心改动】将位置误差乘以 1000 转换为毫米（mm）后，再存入历史矩阵
        err_history(iter, :) = [(e_pos * 1000)', e_ori'];
        
        error_norm = norm([e_pos; e_ori]);
        if error_norm < tol
            tool = true;
            err_history = err_history(1:iter, :); % 截断未使用的行
            return;
        end
        
        % 数值更新 (依旧基于米和弧度更新关节角)
        J = Compute_Jacobian(q, MDH_params);
        lambda = lambda_max * (error_norm + 0.01); 
        A = J' * J + (lambda^2) * eye(7);
        b = J' * [e_pos; e_ori];
        dq = A \ b;
        
        q = q + dq;
    end
end