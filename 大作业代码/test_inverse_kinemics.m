clc; clear; close all;
config_input = [0.2; -0.4; 0.3; 0.6; -0.2; 0.5; 0.1];
T_forward = my_getTransform(config_input);

% === 调试：绘制 psi 在 0~2pi 范围内的残差曲线 ===
psi_test = linspace(0, 2*pi, 360);
res_test = zeros(size(psi_test));

% 这里需要提取出函数内部需要的固定中间变量，以便单独测试子函数
% 临时复制主函数前期的几何提取逻辑：
MDH_params = [0,pi,-0.15643-0.12838;
              0,pi/2,-0.005375-0.006375;
              0,-pi/2,-0.21038-0.21038;
              0,pi/2,-0.006375-0.006375;
              0,-pi/2,-0.20843-0.10593;
              0,pi/2,0;
              0,-pi/2,-0.061525-0.10593;
              0,pi,0];%%EE
d = MDH_params(:,3);
T_7_EE = [1,0,0,0; 0,-1,0,0; 0,0,-1,0; 0,0,0,1];
T_0_7 = T_7_EE/T_forward;
p_w = T_0_7*[0;0;-d(7);1]; P_w = p_w(1:3); P_s = [0;0;-d(1)];
L1 = abs(d(3)); L2 = abs(d(5)); D = norm(P_w-P_s);
cos_val = (L1^2 + L2^2 - D^2) / (2 * L1 * L2);
if cos_val > 1; cos_val = 1; elseif cos_val < -1; cos_val = -1; end
theta4_init = pi - acos(cos_val);
k = (P_w-P_s)/norm(P_w-P_s); Z0 = [0;0;1];
if abs(dot(k,Z0))>0.995
    n_ref = cross(k,[1;0;0])/norm(cross(k,[1;0;0]));
else
    n_ref = cross(k,Z0)/norm(cross(k,Z0));
end

% 循环计算残差
for i = 1:length(psi_test)
    [res_test(i), ~, ~] = resi(psi_test(i), P_w, P_s, MDH_params, theta4_init, n_ref, k);
end

figure;
plot(psi_test, res_test, 'LineWidth', 2);
grid on;
xlabel('\psi (rad)');
ylabel('Residual Error');
title('残差随 \psi 的变化曲线');