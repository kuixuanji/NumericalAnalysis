%% 参数不准确情况下的鲁棒轨迹跟踪控制
% 本脚本复用原大作业中的七自由度机械臂模型和惯性参数。
% 被控对象使用真实参数进行仿真，控制器前馈项故意使用扰动后的
% 质量、质心和惯性张量参数。本模型不计摩擦。
clear; clc; close all;

this_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(this_dir);
addpath(this_dir);
addpath(code_dir);

inertial_true = load_codegen_inertial(fullfile(code_dir, 'MDH_Inertial_Codegen.mat'));
[inertial_hat, uncertainty_info] = make_uncertain_inertial(inertial_true);

ctrl.Kp = diag([120, 115, 100, 85, 55, 35, 25]);
ctrl.Kd = diag([32, 30, 26, 22, 15, 10, 8]);
ctrl.Kr = diag([10, 10, 8, 8, 5, 3, 2]);
ctrl.Lambda = diag([5.0, 5.0, 4.5, 4.5, 3.5, 3.0, 2.5]);
ctrl.phi = 0.08 * ones(7, 1);

sim.T = 8.0;
sim.tspan = linspace(0, sim.T, 401);
sim.options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);

[q_ref0, dq_ref0] = desired_joint_trajectory(0);
q0 = q_ref0 + [0.08; -0.06; 0.05; -0.05; 0.04; -0.03; 0.02];
dq0 = dq_ref0 + zeros(7, 1);
x0 = [q0; dq0];

fprintf('名义模型使用的参数扰动比例：\n');
disp(uncertainty_info);

ode_rhs = @(t, x) closed_loop_rhs(t, x, inertial_true, inertial_hat, ctrl);
[t, x] = ode45(ode_rhs, sim.tspan, x0, sim.options);

n_step = numel(t);
q = x(:, 1:7);
dq = x(:, 8:14);
q_ref = zeros(n_step, 7);
dq_ref = zeros(n_step, 7);
ddq_ref = zeros(n_step, 7);
tau = zeros(n_step, 7);
tau_ff = zeros(n_step, 7);
tau_pd = zeros(n_step, 7);
tau_smc = zeros(n_step, 7);
err_norm = zeros(n_step, 1);

for k = 1:n_step
    [qr, dqr, ddqr] = desired_joint_trajectory(t(k));
    q_ref(k, :) = qr.';
    dq_ref(k, :) = dqr.';
    ddq_ref(k, :) = ddqr.';

    [tau_k, parts] = robust_pd_smc_controller(t(k), q(k, :).', dq(k, :).', ...
        inertial_hat, ctrl);
    tau(k, :) = tau_k.';
    tau_ff(k, :) = parts.tau_ff.';
    tau_pd(k, :) = parts.tau_pd.';
    tau_smc(k, :) = parts.tau_smc.';
    err_norm(k) = norm(qr - q(k, :).');
end

fprintf('最终关节跟踪误差范数：%.6f rad\n', err_norm(end));
fprintf('最大关节跟踪误差范数：%.6f rad\n', max(err_norm));
fprintf('关节跟踪误差均方根范数：%.6f rad\n', sqrt(mean(err_norm.^2)));

plot_tracking_result(t, q, dq, q_ref, dq_ref, tau, tau_ff, tau_pd, tau_smc, err_norm);

function dx = closed_loop_rhs(t, x, inertial_true, inertial_hat, ctrl)
    q = x(1:7);
    dq = x(8:14);

        tau = robust_pd_smc_controller(t, q, dq, inertial_hat, ctrl);

    [M, C, G] = gen3_lagrange_dynamics_param(q, dq, zeros(7, 1), inertial_true);
    ddq = M \ (tau - C * dq - G);

    dx = [dq; ddq];
end

function plot_tracking_result(t, q, dq, q_ref, dq_ref, tau, tau_ff, tau_pd, tau_smc, err_norm)
    figure('Name', '参数不准确情况下的鲁棒跟踪效果', 'Color', 'w');

    subplot(4, 1, 1);
    plot(t, q_ref, '--', 'LineWidth', 1.0); hold on;
    plot(t, q, 'LineWidth', 1.2);
    grid on;
    ylabel('关节角 q / rad');
    title('关节位置跟踪');

    subplot(4, 1, 2);
    plot(t, dq_ref, '--', 'LineWidth', 1.0); hold on;
    plot(t, dq, 'LineWidth', 1.2);
    grid on;
    ylabel('关节角速度 dq / rad/s');
    title('关节速度跟踪');

    subplot(4, 1, 3);
    plot(t, q_ref - q, 'LineWidth', 1.1);
    grid on;
    ylabel('位置误差 e / rad');
    title('关节位置误差');

    subplot(4, 1, 4);
    plot(t, err_norm, 'k', 'LineWidth', 1.5);
    grid on;
    xlabel('时间 t / s');
    ylabel('误差范数 ||e|| / rad');
    title('误差范数');

    figure('Name', '控制力矩组成', 'Color', 'w');
    subplot(4, 1, 1);
    plot(t, tau, 'LineWidth', 1.1);
    grid on;
    ylabel('总力矩 tau');
    title('总控制力矩');

    subplot(4, 1, 2);
    plot(t, tau_ff, 'LineWidth', 1.1);
    grid on;
    ylabel('前馈力矩 tau_{ff}');
    title('名义模型前馈力矩');

    subplot(4, 1, 3);
    plot(t, tau_pd, 'LineWidth', 1.1);
    grid on;
    ylabel('PD力矩 tau_{pd}');
    title('PD反馈力矩');

    subplot(4, 1, 4);
    plot(t, tau_smc, 'LineWidth', 1.1);
    grid on;
    xlabel('时间 t / s');
    ylabel('鲁棒力矩 tau_{smc}');
    title('滑模鲁棒力矩');
end
