%% Lagrange_Dynamics_Main 接口测试脚本
% 本脚本先调用逆运动学求解 q，再将 q、dq、ddq 输入动力学接口。

clear; clc;

this_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(this_dir);
addpath(this_dir);
addpath(code_dir);

MDH_params = [0,  pi,   -0.15643 - 0.12838;
              0,  pi/2, -0.005375 - 0.006375;
              0, -pi/2, -0.21038 - 0.21038;
              0,  pi/2, -0.006375 - 0.006375;
              0, -pi/2, -0.20843 - 0.10593;
              0,  pi/2,  0;
              0, -pi/2, -0.061525 - 0.10593;
              0,  pi,    0];

q_seed = [0.1; -0.2; 0.3; -0.4; 0.2; -0.1; 0.15];
T_target = Forward_Kinematics(q_seed);
q_init = q_seed + 0.02 * ones(7, 1);
[q, ik_success] = Inverse_Kinematics_Numerical(T_target, q_init, MDH_params);

if ~ik_success
    error('Inverse_Kinematics_Numerical 未收敛，测试终止。');
end

qd = [0.2; -0.1; 0.15; 0.05; -0.08; 0.12; -0.04];
qdd = [0.5; -0.3; 0.2; -0.1; 0.25; -0.15; 0.1];

disp('===== 逆运动学 + 动力学接口调用测试 =====');
[M, C, G, tau, T_end] = Lagrange_Dynamics_Main(q, qd, qdd);

check_size(M, [7, 7], 'M');
check_size(C, [7, 7], 'C');
check_size(G, [7, 1], 'G');
check_size(tau, [7, 1], 'tau');
check_size(T_end, [4, 4], 'T_end');
check_finite(M, 'M');
check_finite(C, 'C');
check_finite(G, 'G');
check_finite(tau, 'tau');
check_finite(T_end, 'T_end');

mass_matrix_symmetry_error = max(max(abs(M - M.')));
fprintf('惯性矩阵对称性误差：%.6e\n', mass_matrix_symmetry_error);

if mass_matrix_symmetry_error > 1e-8
    error('惯性矩阵 M 不满足对称性要求。');
end

fprintf('逆运动学末端位置误差：%.6e m\n', norm(T_target(1:3, 4) - T_end(1:3, 4)));
disp('逆运动学 + 动力学接口调用测试通过。');

function check_size(x, expected_size, name)
    if ~isequal(size(x), expected_size)
        error('%s 尺寸错误，应为 %dx%d，实际为 %dx%d。', ...
            name, expected_size(1), expected_size(2), size(x, 1), size(x, 2));
    end
end

function check_finite(x, name)
    if any(~isfinite(x(:)))
        error('%s 中存在 NaN 或 Inf。', name);
    end
end
