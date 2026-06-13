function [M, C, G, tau, T_end] = gen3_lagrange_dynamics_param(q, dq, ddq, inertial)
% 使用外部传入惯性参数的七自由度机械臂拉格朗日动力学。
% 力矩计算式：tau = M(q) * ddq + C(q,dq) * dq + G(q)。

    if nargin < 1 || isempty(q)
        q = zeros(7, 1);
    end
    if nargin < 2 || isempty(dq)
        dq = zeros(7, 1);
    end
    if nargin < 3 || isempty(ddq)
        ddq = zeros(7, 1);
    end
    q = q(:);
    dq = dq(:);
    ddq = ddq(:);
    check_joint_vector(q, 'q');
    check_joint_vector(dq, 'dq');
    check_joint_vector(ddq, 'ddq');

    MDH_params = [0,  pi,   -0.15643 - 0.12838;
                  0,  pi/2, -0.005375 - 0.006375;
                  0, -pi/2, -0.21038 - 0.21038;
                  0,  pi/2, -0.006375 - 0.006375;
                  0, -pi/2, -0.20843 - 0.10593;
                  0,  pi/2,  0;
                  0, -pi/2, -0.061525 - 0.10593;
                  0,  pi,    0];

    g = [0; 0; -9.81];

    M = compute_mass_matrix(q, MDH_params, inertial);
    C = compute_coriolis_matrix(q, dq, MDH_params, inertial);
    G = compute_gravity_vector(q, MDH_params, inertial, g);
    tau = M * ddq + C * dq + G;

    if nargout > 4
        T_end = Forward_Kinematics(q);
    end
end

function M = compute_mass_matrix(q, MDH_params, inertial)
    n = 7;
    M = zeros(n, n);
    [origin, axis_z, T_link] = compute_joint_frames(q, MDH_params);

    for link_id = 1:n
        [m, r_c, I_c] = get_link_inertial(inertial, link_id);
        R_link = T_link{link_id}(1:3, 1:3);
        p_com = T_link{link_id} * [r_c; 1];
        p_com = p_com(1:3);

        Jv = zeros(3, n);
        Jw = zeros(3, n);
        for joint_id = 1:link_id
            z = axis_z(:, joint_id);
            p = origin(:, joint_id);
            Jv(:, joint_id) = cross(z, p_com - p);
            Jw(:, joint_id) = z;
        end

        I_world = R_link * I_c * R_link.';
        M = M + m * (Jv.' * Jv) + Jw.' * I_world * Jw;
    end

    M = 0.5 * (M + M.');
end

function C = compute_coriolis_matrix(q, dq, MDH_params, inertial)
    n = 7;
    h = 1e-6;
    dM = zeros(n, n, n);

    for k = 1:n
        q_plus = q;
        q_minus = q;
        q_plus(k) = q_plus(k) + h;
        q_minus(k) = q_minus(k) - h;

        M_plus = compute_mass_matrix(q_plus, MDH_params, inertial);
        M_minus = compute_mass_matrix(q_minus, MDH_params, inertial);
        dM(:, :, k) = (M_plus - M_minus) / (2 * h);
    end

    C = zeros(n, n);
    for i = 1:n
        for j = 1:n
            c = 0;
            for k = 1:n
                christoffel = 0.5 * (dM(i, j, k) + dM(i, k, j) - dM(j, k, i));
                c = c + christoffel * dq(k);
            end
            C(i, j) = c;
        end
    end
end

function G = compute_gravity_vector(q, MDH_params, inertial, g)
    n = 7;
    h = 1e-6;
    G = zeros(n, 1);

    for k = 1:n
        q_plus = q;
        q_minus = q;
        q_plus(k) = q_plus(k) + h;
        q_minus(k) = q_minus(k) - h;

        V_plus = compute_potential_energy(q_plus, MDH_params, inertial, g);
        V_minus = compute_potential_energy(q_minus, MDH_params, inertial, g);
        G(k) = (V_plus - V_minus) / (2 * h);
    end
end

function V = compute_potential_energy(q, MDH_params, inertial, g)
    n = 7;
    V = 0;
    [~, ~, T_link] = compute_joint_frames(q, MDH_params);

    for link_id = 1:n
        [m, r_c, ~] = get_link_inertial(inertial, link_id);
        p_com = T_link{link_id} * [r_c; 1];
        p_com = p_com(1:3);
        V = V - m * g.' * p_com;
    end
end

function [m, r_c, I_c] = get_link_inertial(inertial, link_id)
    link_data = inertial{link_id};
    m = link_data.m;
    r_c = link_data.r_c(:);
    I_c = link_data.I_c;
end

function [origin, axis_z, T_link] = compute_joint_frames(q, MDH_params)
    n = 7;
    origin = zeros(3, n);
    axis_z = zeros(3, n);
    T_link = cell(1, n);

    T = eye(4);
    for joint_id = 1:n
        a = MDH_params(joint_id, 1);
        alpha = MDH_params(joint_id, 2);
        d = MDH_params(joint_id, 3);
        theta = q(joint_id);

        T = T * Get_MDH_Trans(a, alpha, d, theta);
        T_link{joint_id} = T;
        origin(:, joint_id) = T(1:3, 4);
        axis_z(:, joint_id) = T(1:3, 3);
    end
end

function check_joint_vector(x, name)
    if numel(x) ~= 7
        error('%s 必须包含 7 个关节量。', name);
    end
    if any(~isfinite(x(:)))
        error('%s 中包含 NaN 或 Inf。', name);
    end
end
