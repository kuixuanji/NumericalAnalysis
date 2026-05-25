clear;clc;close all;
sim_result = sim('GEN3_URDF_V12_2024a.slx');
%提取仿真结果
tau_sim = sim_result.torque.signals.values;
t_sim = sim_result.tout;

%调用自己函数
q = zeros(7,length(t_sim));
dq = zeros(7,length(t_sim));
ddq = zeros(7,length(t_sim));
tau_my = zeros(7,length(t_sim));
for i = 1:length(t_sim)
data = load('MDH_Inertial.mat');
MDH_Inertial = data.MDH_Inertial;
[q(:,i), dq(:,i), ddq(:,i)] = trajectory(t_sim(i));
q_current   = q(:, i);
dq_current  = dq(:, i);
ddq_current = ddq(:, i);
[~, ~, ~, tau_my(:,i), ~]=Lagrange_Dynamics_Main(q_current,dq_current,ddq_current);

end

%% 绘图
for i = 1:7
figure;
subplot(2,1,1);
plot(t_sim, tau_sim(:, i), 'b-', 'LineWidth', 2); hold on;
plot(t_sim, tau_my(i,:), 'r--', 'LineWidth', 2); 
grid on;
xlabel('Time (s)', 'FontSize', 11);
ylabel('Torque (N·m)', 'FontSize', 11);
title(['关节',num2str(i), '动力学验证：Simscape vs RNEA'], 'FontSize', 13);
legend('Simscape 物理仿真', 'RNEA 理论计算');
e_i = abs(tau_sim(:, i)-tau_my(i,:)');
subplot(2,1,2);
plot(t_sim,e_i, 'b-', 'LineWidth', 2);
grid on;
xlabel('Time (s)', 'FontSize', 11);
ylabel('Torque Error(N·m)', 'FontSize', 11);
title(['关节',num2str(i), '力矩误差分析'], 'FontSize', 13);
end
function [q, dq, ddq] = trajectory(t)
    % t 由 Clock 模块输入
    A = 0.2;            % 振幅
    w = 2 * pi * 0.3;   % 角频率
    
    % 计算 7x1 的向量
    q   =  A * sin(w * t) * ones(7, 1);       % 位置
    dq  =  A * w * cos(w * t) * ones(7, 1);   % 速度 (q的一阶导)
    ddq = -A * w^2 * sin(w * t) * ones(7, 1); % 加速度 (q的二阶导)
end