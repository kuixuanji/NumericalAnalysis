function result = Forward_Dynamics_Solve(tau_input, q0, dq0, tspan, options)
    if nargin < 4
        error('Usage: result = Forward_Dynamics_Solve(tau_input, q0, dq0, tspan)');
    end
    if nargin < 5 || isempty(options)
        options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
    end

    this_dir = fileparts(mfilename('fullpath'));
    code_dir = fileparts(this_dir);
    addpath(this_dir);
    addpath(code_dir);

    q0 = q0(:);
    dq0 = dq0(:);
    check_joint_vector(q0, 'q0');
    check_joint_vector(dq0, 'dq0');
    validate_tspan(tspan);

    tau_fun = make_torque_function(tau_input, tspan);
    x0 = [q0; dq0];

    ode_rhs = @(t, x) forward_dynamics_rhs(t, x, tau_fun);
    [t, x] = ode45(ode_rhs, tspan(:), x0, options);

    n_step = numel(t);
    q = x(:, 1:7);
    dq = x(:, 8:14);
    ddq = zeros(n_step, 7);
    tau = zeros(n_step, 7);
    T_end = zeros(4, 4, n_step);

    for k = 1:n_step
        qk = q(k, :).';
        dqk = dq(k, :).';
        tauk = tau_fun(t(k), qk, dqk);
        ddq(k, :) = compute_forward_acceleration(qk, dqk, tauk).';
        tau(k, :) = tauk.';
        T_end(:, :, k) = Forward_Kinematics(qk);
    end

    result.t = t;
    result.q = q;
    result.dq = dq;
    result.ddq = ddq;
    result.tau = tau;
    result.T_end = T_end;
end

function dx = forward_dynamics_rhs(t, x, tau_fun)
    q = x(1:7);
    dq = x(8:14);
    tau = tau_fun(t, q, dq);
    ddq = compute_forward_acceleration(q, dq, tau);
    dx = [dq; ddq];
end

function ddq = compute_forward_acceleration(q, dq, tau)
    [M, C, G] = Lagrange_Dynamics_Main(q, dq, zeros(7, 1));
    ddq = M \ (tau(:) - C * dq(:) - G);
end

function tau_fun = make_torque_function(tau_input, tspan)
    if isa(tau_input, 'function_handle')
        tau_fun = @(t, q, dq) validate_tau(tau_input(t, q, dq));
        return;
    end

    if isnumeric(tau_input) && isequal(size(tau_input), [7, 1])
        tau_const = tau_input(:);
        tau_fun = @(~, ~, ~) tau_const;
        return;
    end

    if isnumeric(tau_input) && isequal(size(tau_input), [1, 7])
        tau_const = tau_input(:);
        tau_fun = @(~, ~, ~) tau_const;
        return;
    end

    if isnumeric(tau_input) && size(tau_input, 2) == 7 && size(tau_input, 1) == numel(tspan)
        t_samples = tspan(:);
        tau_samples = tau_input;
        tau_fun = @(t, ~, ~) interp1(t_samples, tau_samples, t, 'linear', 'extrap').';
        return;
    end

    error(['tau_input must be a 7x1 vector, a 1x7 vector, ', ...
           'a numel(tspan)x7 matrix, or a function handle.']);
end

function tau = validate_tau(tau)
    tau = tau(:);
    check_joint_vector(tau, 'tau');
end

function check_joint_vector(x, name)
    if numel(x) ~= 7
        error('%s must include 7 joints', name);
    end
    if any(~isfinite(x(:)))
        error('%s contains NaN or Inf', name);
    end
end

function validate_tspan(tspan)
    if ~isnumeric(tspan) || numel(tspan) < 2
        error('tspan must contain at least two time points');
    end
    if any(~isfinite(tspan(:))) || any(diff(tspan(:)) <= 0)
        error('tspan must be finite and strictly increasing');
    end
end
