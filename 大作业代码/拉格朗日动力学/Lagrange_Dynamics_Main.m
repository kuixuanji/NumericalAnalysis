function [M, C, G, tau, T_end] = Lagrange_Dynamics_Main(q, qd, qdd)
% 调用格式：由已知关节状态直接计算动力学
% [M, C, G, tau, T_end] = Lagrange_Dynamics_Main(q, qd, qdd)
%
% 输入：
% q        7x1 关节角（由外部逆运动学求解器预先给定）
% qd       7x1 关节角速度
% qdd      7x1 关节角加速度
%
% 输出：
% M          7x7 惯性矩阵
% C          7x7 科氏力和离心力矩阵
% G          7x1 重力项
% tau        7x1 关节驱动力矩
% T_end      4x4 当前末端位姿

    if nargin < 1 || isempty(q)
        q = zeros(7, 1);
    end
    if nargin < 2 || isempty(qd)
        qd = zeros(7, 1);
    end
    if nargin < 3 || isempty(qdd)
        qdd = zeros(7, 1);
    end

    q = q(:);
    qd = qd(:);
    qdd = qdd(:);
    check_joint_vector(q, 'q');
    check_joint_vector(qd, 'qd');
    check_joint_vector(qdd, 'qdd');

    this_dir = fileparts(mfilename('fullpath'));
    code_dir = fileparts(this_dir);
    addpath(code_dir);

    load(fullfile(code_dir, 'MDH_Inertial.mat'), 'MDH_Inertial');

    MDH_params = [0,  pi,   -0.15643 - 0.12838;
                  0,  pi/2, -0.005375 - 0.006375;
                  0, -pi/2, -0.21038 - 0.21038;
                  0,  pi/2, -0.006375 - 0.006375;
                  0, -pi/2, -0.20843 - 0.10593;
                  0,  pi/2,  0;
                  0, -pi/2, -0.061525 - 0.10593;
                  0,  pi,    0];      % 末端执行器固定变换

    g = [0; 0; -9.81];
    % 正运动学按要求调用外部 Forward_Kinematics.m。
    T_end = Forward_Kinematics(q);

    % 拉格朗日动力学模型：tau = M(q) * qdd + C(q, qd) * qd + G(q)。
    [M, C, G, tau] = compute_lagrange_dynamics(q, qd, qdd, MDH_params, MDH_Inertial, g);
end

function check_joint_vector(x, name)
    if numel(x) ~= 7
        error('%s must include 7 joints', name);
    end
end

function [M, C, G, tau] = compute_lagrange_dynamics(q, qd, qdd, MDH_params, inertial, g)
    q = q(:);
    qd = qd(:);
    qdd = qdd(:);

    M = compute_mass_matrix(q, MDH_params, inertial);
    C = compute_coriolis_matrix(q, qd, MDH_params, inertial);
    G = compute_gravity_vector(q, MDH_params, inertial, g);

    tau = M * qdd + C * qd + G;
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

function C = compute_coriolis_matrix(q, qd, MDH_params, inertial)
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
                c = c + christoffel * qd(k);
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
