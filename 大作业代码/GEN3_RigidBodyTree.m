%% Kinova Gen3 7-DOF 路径验证 (自动匹配节点名)
clear; clc;

% 1. 加载模型
try
    robot = importrobot('GEN3_URDF_V12');
catch
    error('未找到 GEN3_URDF_V12 文件，请确保它在当前文件夹或路径中。');
end
robot.DataFormat = 'column';

% 2. 自动定位末端执行器名称
% 通常最后一个 body 就是末端，或者你可以根据 robot.BodyNames 手动指定
all_bodies = robot.BodyNames;
eeName = all_bodies{end}; 
fprintf('正在对比基准模型末端: [%s]\n', eeName);

% 3. 生成测试路径
num_steps = 50;
q_start = zeros(7, 1);
q_end = [3; 0.4; -0.1; 3; 0; 0.7; 3]; % 随机一个测试姿态

t = linspace(0, 1, num_steps);
q_path = q_start*(1-t) + q_end*t;

pos_custom = zeros(3, num_steps);
pos_robot  = zeros(3, num_steps);

% 4. 计算对比
fprintf('开始验证...\n');
for i = 1:num_steps
    q = q_path(:, i);
    
    % 自定义函数
    T_c = my_getTransform(q);
    pos_custom(:, i) = T_c(1:3, 4);
    
    % 官方函数
    T_r = getTransform(robot, q, eeName);
    pos_robot(:, i) = T_r(1:3, 4);
end

% 5. 绘图与误差分析
err = sqrt(sum((pos_custom - pos_robot).^2, 1)) * 1000; % 误差(mm)

figure('Color', 'w');
subplot(1,2,1);
plot3(pos_custom(1,:), pos_custom(2,:), pos_custom(3,:), 'r', 'LineWidth', 3); hold on;
plot3(pos_robot(1,:), pos_robot(2,:), pos_robot(3,:), 'c--', 'LineWidth', 2);
grid on; axis equal; view(3);
title('轨迹对比'); legend('Custom MDH', 'URDF Model');

subplot(1,2,2);
plot(err, 'r', 'LineWidth', 1.5);
title('位置误差 (mm)'); ylabel('Error (mm)'); xlabel('Step');
grid on;

fprintf('最大位置偏差: %.6f mm\n', max(err));
