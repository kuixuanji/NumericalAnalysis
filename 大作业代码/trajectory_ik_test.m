clc; clear; close all;
%% 1. 配置机械臂的 MDH 参数
MDH_params = [
            0,      pi,      -0.1564-0.1284;  
            0,      pi/2,    -0.0054-0.0064;  
            0,     -pi/2,    -0.2104-0.2104;  
            0,      pi/2,    -0.0064-0.0064;  
            0,     -pi/2,    -0.2084-0.1059;  
            0,      pi/2,     0;              
            0,     -pi/2,    -0.0615-0.1059;  
            0,      pi,       0               
            ];
%% 2. 生成笛卡尔空间运动轨迹 (以 YZ 平面的圆为例)
num_points = 200; % 轨迹总插补点数
theta = linspace(0, 2*pi, num_points);

% 定义圆心和半径 (单位: 米)
center_x = 0.48; 
center_y = 0.0;
center_z = 0.25;
radius   = 0.3;
    
% 期望的末端姿态
R_target = [1,  0,  0;
            0,  -1,  0;
            0,  0,  1];
            
% 预分配矩阵存储结果
p_actual_history = zeros(num_points, 3);
q_history = zeros(num_points, 7);
error_history = zeros(num_points, 1); 

%% 3. 连续逆运动学求解 (核心：使用 Warm Start)
% 初始点盲解猜想
q_current = zeros(7, 1); 

fprintf('开始追踪连续空间运动曲线...\n');
tic;
for i = 1:num_points
    p_target = [center_x ; 
                center_y + radius * cos(theta(i)); 
                center_z + radius * sin(theta(i))];
            
    T_target = eye(4);
    T_target(1:3, 1:3) = R_target;
    T_target(1:3, 4)   = p_target;
    
    % 正常的连续追踪
    [q_sol, ~] = Inverse_Kinematics_Numerical(T_target, q_current, MDH_params);
        
    q_current = q_sol; 
    q_history(i, :) = q_sol';
        
    % 正运动学计算实际达到的位姿
    T_act = Forward_Kinematics(q_sol);
    p_actual = T_act(1:3, 4)';
    p_actual_history(i, :) = p_actual;
        
    % 计算期望轨迹与实际轨迹之间的欧氏距离误差 (单位：毫米 mm)
    error_history(i) = norm(p_target' - p_actual) * 1000; 
end
toc;
    
%% 4. 绘制测试与评估曲线 
figure('Name', '机械臂连续轨迹逆运动学测试与精度评估', 'NumberTitle', 'off', 'Position', [50, 100, 1500, 450]);

% ---- 子图 1：笛卡尔空间三维轨迹追踪效果 ----
subplot(1, 3, 1);
p_desired = [repmat(center_x, 1, num_points); center_y + radius * cos(theta); center_z + radius * sin(theta)]';
plot3(p_desired(:,1), p_desired(:,2), p_desired(:,3), 'r--', 'LineWidth', 2, 'DisplayName', 'Desired Path'); hold on;
plot3(p_actual_history(:,1), p_actual_history(:,2), p_actual_history(:,3), 'b-', 'LineWidth', 1.5, 'DisplayName', 'IK Tracked Path');
grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('末端三维轨迹追踪对比');
legend('Location', 'best');
    
% ---- 子图 2：关节空间连续运动曲线 ----
subplot(1, 3, 2);
plot(1:num_points, q_history, 'LineWidth', 1.5);
grid on;
xlabel('Trajectory Points'); ylabel('Joint Angles (rad)');
title('7个关节的连续运动角曲线');
legend('q1','q2','q3','q4','q5','q6','q7', 'Location', 'best');
    
% ---- 子图 3：追踪误差距离曲线 ----
subplot(1, 3, 3);
plot(1:num_points, error_history, 'm-', 'LineWidth', 2);
grid on;
xlabel('Trajectory Points'); ylabel('Position Error (mm)');
title('末端轨迹位置误差 (Tracking Error)');
    
% 在图上标注最大误差和平均误差
max_err = max(error_history);
mean_err = mean(error_history);
txt = sprintf('Max Error: %.4f mm\nMean Error: %.4f mm', max_err, mean_err);
annotation('textbox', [0.76, 0.75, 0.13, 0.12], 'String', txt, 'FitBoxToText', 'on', 'BackgroundColor', 'w');

