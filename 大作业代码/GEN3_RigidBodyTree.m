%% Kinova Gen3 7-DOF 路径验证（自动匹配 URDF）
clear; clc;

% 1) 加载模型：优先加载脚本同目录下的 GEN3_URDF_V12.urdf
script_dir = fileparts(mfilename('fullpath'));
candidate_urdf = fullfile(script_dir, 'GEN3_URDF_V12.urdf');

if exist(candidate_urdf, 'file')
    urdf_path = candidate_urdf;
else
    % 2) 兜底：递归搜索同目录及子目录下的 urdf
    urdf_list = dir(fullfile(script_dir, '**', '*.urdf'));
    if isempty(urdf_list)
        error('未找到 URDF 文件。请确认 GEN3_URDF_V12.urdf 在当前目录或子目录中。');
    end

    % 优先匹配包含 GEN3_URDF_V12 的文件名
    idx = find(contains({urdf_list.name}, 'GEN3_URDF_V12', 'IgnoreCase', true), 1);
    if isempty(idx)
        idx = 1;
    end
    urdf_path = fullfile(urdf_list(idx).folder, urdf_list(idx).name);
    fprintf('未直接找到 GEN3_URDF_V12.urdf，改用：%s\n', urdf_path);
end

robot = importrobot(urdf_path);
robot.DataFormat = 'column';

% 2. 自动定位末端执行器名称
all_bodies = robot.BodyNames;
eeName = all_bodies{end};
fprintf('正在对比基准模型末端: [%s]\n', eeName);

% 3. 生成测试路径
num_steps = 50;
q_start = zeros(7, 1);
q_end = [3; 0.4; -0.1; 3; 0; 0.7; 3];

% 线性插值每个关节角
t = linspace(0, 1, num_steps);
q_path = q_start*(1-t) + q_end*t;

pos_custom = zeros(3, num_steps);
pos_robot  = zeros(3, num_steps);

% 4. 计算对比
fprintf('开始验证...\n');
for i = 1:num_steps
    q = q_path(:, i);

    % 自定义函数（如果你项目里用的是 Forward_Kinematics，可替换这里）
    %T_c = my_getTransform(q);
    T_c = Forward_Kinematics(q);
    pos_custom(:, i) = T_c(1:3, 4);

    % 官方函数
    T_r = getTransform(robot, q, eeName);
    pos_robot(:, i) = T_r(1:3, 4);
end

% 5. 绘图与误差分析
err = sqrt(sum((pos_custom - pos_robot).^2, 1)) * 1000; % 误差 mm

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
