function [tau, parts] = robust_pd_smc_controller(t, q, dq, inertial_hat, ctrl)
% 前馈 + PD + 饱和滑模鲁棒项控制器。

    [q_ref, dq_ref, ddq_ref] = desired_joint_trajectory(t);

    [~, ~, ~, tau_model] = gen3_lagrange_dynamics_param( ...
        q_ref, dq_ref, ddq_ref, inertial_hat);

    e = q_ref - q(:);
    de = dq_ref - dq(:);
    s = de + ctrl.Lambda * e;

    tau_ff = tau_model;
    tau_pd = ctrl.Kp * e + ctrl.Kd * de;
    tau_smc = ctrl.Kr * sat_vector(s ./ ctrl.phi);
    tau = tau_ff + tau_pd + tau_smc;

    if nargout > 1
        parts.q_ref = q_ref;
        parts.dq_ref = dq_ref;
        parts.ddq_ref = ddq_ref;
        parts.e = e;
        parts.de = de;
        parts.s = s;
        parts.tau_ff = tau_ff;
        parts.tau_pd = tau_pd;
        parts.tau_smc = tau_smc;
    end
end
