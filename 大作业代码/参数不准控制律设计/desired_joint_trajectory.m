function [q_ref, dq_ref, ddq_ref] = desired_joint_trajectory(t)
% 平滑有界的关节空间期望轨迹。

    q_offset = [0.00; -0.35; 0.20; -0.70; 0.15; 0.45; 0.05];
    amp = [0.22; 0.18; 0.16; 0.14; 0.12; 0.10; 0.08];
    omega = [0.70; 0.60; 0.80; 0.55; 0.90; 0.75; 0.65];
    phase = [0.00; 0.40; 0.80; 1.20; 0.60; 1.00; 1.40];

    q_ref = q_offset + amp .* sin(omega * t + phase);
    dq_ref = amp .* omega .* cos(omega * t + phase);
    ddq_ref = -amp .* omega.^2 .* sin(omega * t + phase);
end
