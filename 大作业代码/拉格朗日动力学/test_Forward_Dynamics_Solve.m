%% Forward_Dynamics_Solve interface test script
% Given tau, q0 and dq0, integrate the Lagrange forward dynamics.

clear; clc; close all;

this_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(this_dir);
addpath(this_dir);
addpath(code_dir);

q0 = [0.10; -0.20; 0.30; -0.40; 0.20; -0.10; 0.15];
dq0 = zeros(7, 1);
tspan = linspace(0, 2, 101);

% Example 1: constant torque.
tau_const = [0.8; -0.5; 0.4; -0.2; 0.15; -0.1; 0.05];
result_const = Forward_Dynamics_Solve(tau_const, q0, dq0, tspan);

% Example 2: time-varying torque. Uncomment to use.
% tau_fun = @(t, q, dq) tau_const + 0.1 * sin(2 * pi * t) * ones(7, 1);
% result_fun = Forward_Dynamics_Solve(tau_fun, q0, dq0, tspan);

fprintf('Final q:\n');
disp(result_const.q(end, :).');
fprintf('Final dq:\n');
disp(result_const.dq(end, :).');
fprintf('Final ddq:\n');
disp(result_const.ddq(end, :).');

figure('Name', 'Forward Dynamics Result');
subplot(3, 1, 1);
plot(result_const.t, result_const.q, 'LineWidth', 1.1);
grid on;
ylabel('q / rad');
title('Joint Position');

subplot(3, 1, 2);
plot(result_const.t, result_const.dq, 'LineWidth', 1.1);
grid on;
ylabel('dq / rad/s');
title('Joint Velocity');

subplot(3, 1, 3);
plot(result_const.t, result_const.ddq, 'LineWidth', 1.1);
grid on;
xlabel('t / s');
ylabel('ddq / rad/s^2');
title('Joint Acceleration');
