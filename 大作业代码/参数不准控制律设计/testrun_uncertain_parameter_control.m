%% 参数不准确情况下的鲁棒轨迹跟踪控制：诊断版
% 本脚本复用原大作业中的七自由度机械臂模型和惯性参数。
% 被控对象使用真实参数进行仿真，控制器前馈项故意使用扰动后的
% 质量、质心和惯性张量参数。本模型不计摩擦。
% 诊断版保留原控制律结构，并额外加入力矩限幅、关节限位、状态监测和奇异矩阵检查。
clear; clc; close all;

this_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(this_dir);
addpath(this_dir);
addpath(code_dir);

data = load(fullfile(code_dir, 'MDH_Inertial.mat'), 'MDH_Inertial');
inertial_true = data.MDH_Inertial;
[inertial_hat, uncertainty_info] = make_uncertain_inertial(inertial_true);

test.use_true_feedforward = false;
test.gain_scale = 0.10;
test.initial_error_scale = 0.25;
test.feedback_ramp_time = 0.5;
test.joint2_kp_scale = 1.5;
test.joint2_kd_scale = 3.5;
test.joint2_kr_scale = 0.3;

if test.use_true_feedforward
    inertial_ctrl = inertial_true;
else
    inertial_ctrl = inertial_hat;
end

ctrl.Kp = test.gain_scale * diag([120, 115, 100, 85, 55, 35, 25]);
ctrl.Kd = test.gain_scale * diag([32, 30, 26, 22, 15, 10, 8]);
ctrl.Kr = test.gain_scale * diag([10, 10, 8, 8, 5, 3, 2]);
ctrl.Kp(2, 2) = test.joint2_kp_scale * ctrl.Kp(2, 2);
ctrl.Kd(2, 2) = test.joint2_kd_scale * ctrl.Kd(2, 2);
ctrl.Kr(2, 2) = test.joint2_kr_scale * ctrl.Kr(2, 2);
ctrl.Lambda = diag([5.0, 5.0, 4.5, 4.5, 3.5, 3.0, 2.5]);
ctrl.phi = 0.2 * ones(7, 1);
ctrl.tau_limit = [39; 39; 39; 39; 9; 9; 9];
ctrl.feedback_ramp_time = test.feedback_ramp_time;
ctrl.q_lower = [-inf; -2.24; -inf; -2.57; -inf; -2.09; -inf];
ctrl.q_upper = [ inf;  2.24;  inf;  2.57;  inf;  2.09;  inf];
ctrl.motor_inertia = 0.05 * ones(7, 1);
ctrl.q_limit_margin = 0.35;
ctrl.q_limit_Kp = [0; 30; 0; 15; 0; 8; 0];
ctrl.q_limit_Kd = [0; 8; 0; 4; 0; 2; 0];

sim.T = 8.0;
sim.dt = 0.001;
sim.tspan = 0:sim.dt:sim.T;

[q_ref0, dq_ref0] = desired_joint_trajectory(0);
q0 = q_ref0 + test.initial_error_scale * ...
    [0.08; -0.06; 0.05; -0.05; 0.04; -0.03; 0.02];
dq0 = dq_ref0 + zeros(7, 1);
x0 = [q0; dq0];

fprintf('名义模型使用的参数扰动比例：\n');
disp(uncertainty_info);
if test.use_true_feedforward
    fprintf('诊断设置：前馈项使用真实惯性参数，对照参数扰动影响。\n');
else
    fprintf('诊断设置：前馈项使用扰动后的名义惯性参数。\n');
end
fprintf('诊断设置：gain_scale=%.3f, initial_error_scale=%.3f, feedback_ramp_time=%.3f s\n', ...
    test.gain_scale, test.initial_error_scale, test.feedback_ramp_time);
fprintf('关节2额外调节：Kp x %.3f, Kd x %.3f, Kr x %.3f, phi=%.3f\n', ...
    test.joint2_kp_scale, test.joint2_kd_scale, test.joint2_kr_scale, ctrl.phi(1));
fprintf('前向动力学求解加入等效电机惯量：diag(Jm)=%.3f kg*m^2\n', ctrl.motor_inertia(1));
fprintf('关节限位软保护：margin=%.3f rad\n', ctrl.q_limit_margin);

ode_rhs = @(t, x) closed_loop_rhs(t, x, inertial_true, inertial_ctrl, ctrl);
[t, x, stop_info] = fixed_step_diagnostic_integrate( ...
    ode_rhs, sim.tspan, x0, inertial_true, inertial_ctrl, ctrl);

if stop_info.stopped
    fprintf('\n诊断版固定步长仿真提前终止：t = %.6f s\n', stop_info.t);
    fprintf('终止原因：%s\n', stop_info.reason);
    report_state_diagnostics(stop_info.t, stop_info.x, inertial_true, inertial_ctrl, ctrl);
end

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

    [tau_raw, parts] = diagnostic_controller(t(k), q(k, :).', dq(k, :).', ...
        inertial_ctrl, ctrl);
    tau_k = limit_tau(tau_raw, ctrl.tau_limit);
    tau(k, :) = tau_k.';
    tau_ff(k, :) = parts.tau_ff.';
    tau_pd(k, :) = parts.tau_pd.';
    tau_smc(k, :) = parts.tau_smc.';
    err_norm(k) = norm(qr - q(k, :).');
end

fprintf('最终关节跟踪误差范数：%.6f rad\n', err_norm(end));
fprintf('最大关节跟踪误差范数：%.6f rad\n', max(err_norm));
fprintf('关节跟踪误差均方根范数：%.6f rad\n', sqrt(mean(err_norm.^2)));
lower_violation = max(max(ctrl.q_lower.' - q));
upper_violation = max(max(q - ctrl.q_upper.'));
joint_limit_violation = max([0, lower_violation, upper_violation]);
fprintf('最大关节位置限位超出量：%.6f rad\n', joint_limit_violation);

plot_tracking_result(t, q, dq, q_ref, dq_ref, tau, tau_ff, tau_pd, tau_smc, err_norm);

function dx = closed_loop_rhs(t, x, inertial_true, inertial_hat, ctrl)
    q = x(1:7);
    dq = x(8:14);

    tau_raw = diagnostic_controller(t, q, dq, inertial_hat, ctrl);
    tau = limit_tau(tau_raw, ctrl.tau_limit);

    [M, C, G] = gen3_lagrange_dynamics_param(q, dq, zeros(7, 1), inertial_true);
    rhs = tau - C * dq - G;
    M_solve = regularize_mass_matrix(M, ctrl);
    info = compute_diagnostics(t, q, dq, tau, M, rhs, M_solve);
    assert_diagnostics_ok(info);

    ddq = M_solve \ rhs;
    if any(~isfinite(ddq))
        error('t=%.6f: ddq 中存在 NaN 或 Inf。', t);
    end

    dx = [dq; ddq];
end

function [t_out, x_out, stop_info] = fixed_step_diagnostic_integrate( ...
        ode_rhs, tspan, x0, inertial_true, inertial_hat, ctrl)
    n_step = numel(tspan);
    x_all = zeros(n_step, numel(x0));
    x_all(1, :) = x0(:).';

    stop_info.stopped = false;
    stop_info.t = tspan(end);
    stop_info.x = x0(:);
    stop_info.reason = '';

    fprintf('\n开始固定步长诊断积分，共 %d 步。积分格式：半隐式 Euler。\n', n_step - 1);
    report_state_diagnostics(tspan(1), x0(:), inertial_true, inertial_hat, ctrl);

    for k = 1:n_step - 1
        t = tspan(k);
        h = tspan(k + 1) - tspan(k);
        xk = x_all(k, :).';

        try
            dx = ode_rhs(t, xk);
        catch ME
            stop_info.stopped = true;
            stop_info.t = t;
            stop_info.x = xk;
            stop_info.reason = ME.message;
            t_out = tspan(1:k).';
            x_out = x_all(1:k, :);
            return;
        end

        if any(~isfinite(dx))
            stop_info.stopped = true;
            stop_info.t = t;
            stop_info.x = xk;
            stop_info.reason = 'dx 中存在 NaN 或 Inf。';
            t_out = tspan(1:k).';
            x_out = x_all(1:k, :);
            return;
        end

        q_next = xk(1:7) + h * (xk(8:14) + h * dx(8:14));
        dq_next = xk(8:14) + h * dx(8:14);
        x_next = [q_next; dq_next];
        x_all(k + 1, :) = x_next.';

        [is_bad, reason] = is_bad_state(tspan(k + 1), x_next, ...
            inertial_true, inertial_hat, ctrl);
        if is_bad
            stop_info.stopped = true;
            stop_info.t = tspan(k + 1);
            stop_info.x = x_next;
            stop_info.reason = reason;
            t_out = tspan(1:k + 1).';
            x_out = x_all(1:k + 1, :);
            return;
        end

        if mod(k, 500) == 0 || k == n_step - 1
            fprintf('  step %4d/%4d, t = %.3f s\n', k, n_step - 1, t);
            report_state_diagnostics(t, xk, inertial_true, inertial_hat, ctrl);
        end
    end

    t_out = tspan(:);
    x_out = x_all;
end

function [is_bad, reason] = is_bad_state(t, x, inertial_true, inertial_hat, ctrl)
    is_bad = false;
    reason = '';

    if any(~isfinite(x))
        is_bad = true;
        reason = '状态 x 中存在 NaN 或 Inf。';
        return;
    end

    q = x(1:7);
    dq = x(8:14);

    lower_violation = q < ctrl.q_lower;
    upper_violation = q > ctrl.q_upper;
    if any(lower_violation) || any(upper_violation)
        is_bad = true;
        lower_amount = max(ctrl.q_lower - q);
        upper_amount = max(q - ctrl.q_upper);
        reason = sprintf('关节位置越界，最大超出量 %.3e rad。', ...
            max([lower_amount, upper_amount]));
        return;
    end

    tau_raw = diagnostic_controller(t, q, dq, inertial_hat, ctrl);
    tau = limit_tau(tau_raw, ctrl.tau_limit);
    [M, C, G] = gen3_lagrange_dynamics_param(q, dq, zeros(7, 1), inertial_true);
    rhs = tau - C * dq - G;
    M_solve = regularize_mass_matrix(M, ctrl);
    info = compute_diagnostics(t, q, dq, tau, M, rhs, M_solve);

    if ~info.is_finite
        is_bad = true;
        reason = '状态、力矩或动力学矩阵中存在 NaN/Inf。';
    elseif info.norm_q > 20 || info.norm_dq > 100
        is_bad = true;
        reason = sprintf('状态发散，norm(q)=%.3e, norm(dq)=%.3e。', ...
            info.norm_q, info.norm_dq);
    elseif info.rcond_M < 1e-10 || info.min_eig_M <= 0
        is_bad = true;
        reason = sprintf('M 病态或非正定，rcond(M)=%.3e, minEig(M)=%.3e。', ...
            info.rcond_M, info.min_eig_M);
    end
end

function tau = limit_tau(tau, tau_limit)
    tau = max(min(tau(:), tau_limit(:)), -tau_limit(:));
end

function tau_guard = joint_limit_guard(q, dq, ctrl)
    tau_guard = zeros(7, 1);

    for i = 1:7
        if isfinite(ctrl.q_lower(i))
            lower_zone = ctrl.q_lower(i) + ctrl.q_limit_margin;
            if q(i) < lower_zone
                pos_err = lower_zone - q(i);
                vel_toward_limit = max(-dq(i), 0);
                tau_guard(i) = tau_guard(i) + ...
                    ctrl.q_limit_Kp(i) * pos_err + ...
                    ctrl.q_limit_Kd(i) * vel_toward_limit;
            end
        end

        if isfinite(ctrl.q_upper(i))
            upper_zone = ctrl.q_upper(i) - ctrl.q_limit_margin;
            if q(i) > upper_zone
                pos_err = q(i) - upper_zone;
                vel_toward_limit = max(dq(i), 0);
                tau_guard(i) = tau_guard(i) - ...
                    ctrl.q_limit_Kp(i) * pos_err - ...
                    ctrl.q_limit_Kd(i) * vel_toward_limit;
            end
        end
    end
end

function M_solve = regularize_mass_matrix(M, ctrl)
    M_solve = 0.5 * (M + M.') + diag(ctrl.motor_inertia);
end

function [tau, parts] = diagnostic_controller(t, q, dq, inertial_hat, ctrl)
    [tau_full, parts] = robust_pd_smc_controller(t, q, dq, inertial_hat, ctrl);

    if ctrl.feedback_ramp_time > 0
        feedback_scale = min(max(t / ctrl.feedback_ramp_time, 0), 1);
    else
        feedback_scale = 1;
    end

    parts.tau_pd_raw = parts.tau_pd;
    parts.tau_smc_raw = parts.tau_smc;
    parts.tau_pd = feedback_scale * parts.tau_pd_raw;
    parts.tau_smc = feedback_scale * parts.tau_smc_raw;
    parts.tau_limit_guard = joint_limit_guard(q(:), dq(:), ctrl);
    tau = parts.tau_ff + parts.tau_pd + parts.tau_smc + parts.tau_limit_guard;
    parts.feedback_scale = feedback_scale;
    parts.tau_full = tau_full;
end

function info = compute_diagnostics(t, q, dq, tau, M, rhs, M_solve)
    if nargin < 7
        M_solve = M;
    end

    M_sym = 0.5 * (M + M.');
    M_solve_sym = 0.5 * (M_solve + M_solve.');
    ddq_raw = M \ rhs;
    ddq = M_solve \ rhs;

    info.t = t;
    info.norm_q = norm(q);
    info.norm_dq = norm(dq);
    info.norm_tau = norm(tau);
    info.norm_rhs = norm(rhs);
    info.norm_ddq_raw = norm(ddq_raw);
    info.norm_ddq = norm(ddq);
    info.max_abs_tau = max(abs(tau));
    info.max_abs_rhs = max(abs(rhs));
    info.max_abs_ddq_raw = max(abs(ddq_raw));
    info.max_abs_ddq = max(abs(ddq));
    info.rcond_M = rcond(M_sym);
    info.min_eig_M = min(eig(M_sym));
    info.rcond_M_solve = rcond(M_solve_sym);
    info.min_eig_M_solve = min(eig(M_solve_sym));
    info.is_finite = all(isfinite(q)) && all(isfinite(dq)) && ...
        all(isfinite(tau)) && all(isfinite(M(:))) && ...
        all(isfinite(rhs)) && all(isfinite(ddq_raw)) && all(isfinite(ddq));
end

function assert_diagnostics_ok(info)
    if ~info.is_finite
        error('t=%.6f: 状态、力矩或动力学矩阵中存在 NaN/Inf。', info.t);
    end
    if info.norm_q > 20 || info.norm_dq > 100
        error(['t=%.6f: 状态发散，norm(q)=%.3e, norm(dq)=%.3e, ', ...
            'norm(tau)=%.3e。'], info.t, info.norm_q, info.norm_dq, info.norm_tau);
    end
    if info.rcond_M < 1e-10 || info.min_eig_M <= 0
        error(['t=%.6f: M 病态或非正定，rcond(M)=%.3e, ', ...
            'minEig(M)=%.3e, norm(q)=%.3e, norm(dq)=%.3e, norm(tau)=%.3e。'], ...
            info.t, info.rcond_M, info.min_eig_M, ...
            info.norm_q, info.norm_dq, info.norm_tau);
    end
end

function report_state_diagnostics(t, x, inertial_true, inertial_hat, ctrl)
    q = x(1:7);
    dq = x(8:14);
    [tau_raw, parts] = diagnostic_controller(t, q, dq, inertial_hat, ctrl);
    tau = limit_tau(tau_raw, ctrl.tau_limit);
    [M, C, G] = gen3_lagrange_dynamics_param(q, dq, zeros(7, 1), inertial_true);
    rhs = tau - C * dq - G;
    M_solve = regularize_mass_matrix(M, ctrl);
    info = compute_diagnostics(t, q, dq, tau, M, rhs, M_solve);

    fprintf('  norm(q)    = %.6e\n', info.norm_q);
    fprintf('  norm(dq)   = %.6e\n', info.norm_dq);
    fprintf('  norm(tau raw / limited) = %.6e / %.6e\n', norm(tau_raw), info.norm_tau);
    fprintf('  max |tau raw / limited| = %.6e / %.6e\n', max(abs(tau_raw)), info.max_abs_tau);
    fprintf('  norm(tau_limit_guard) = %.6e, joint2 guard = %.6e\n', ...
        norm(parts.tau_limit_guard), parts.tau_limit_guard(2));
    fprintf('  norm(rhs)  = %.6e\n', info.norm_rhs);
    fprintf('  norm(ddq raw / solve) = %.6e / %.6e\n', info.norm_ddq_raw, info.norm_ddq);
    fprintf('  max |ddq raw / solve| = %.6e / %.6e\n', info.max_abs_ddq_raw, info.max_abs_ddq);
    fprintf('  rcond(M)   = %.6e\n', info.rcond_M);
    fprintf('  minEig(M)  = %.6e\n', info.min_eig_M);
    fprintf('  rcond(M+Jm)= %.6e\n', info.rcond_M_solve);
    fprintf('  minEig(M+Jm)= %.6e\n', info.min_eig_M_solve);
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

    legend;
end
