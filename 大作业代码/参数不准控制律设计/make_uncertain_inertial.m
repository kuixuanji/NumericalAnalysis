function [inertial_hat, info] = make_uncertain_inertial(inertial_true)
% 构造固定偏差的名义惯性参数模型。
% 参数扰动使用固定比例，保证每次实验结果可以复现。

    mass_scale = [0.80; 1.20; 0.85; 1.15; 0.90; 1.25; 0.75];
    com_scale = [1.10; 0.90; 1.08; 0.92; 1.12; 0.88; 1.05];
    inertia_scale = [0.78; 1.18; 0.82; 1.22; 0.88; 1.15; 0.80];

    inertial_hat = inertial_true;
    for i = 1:7
        inertial_hat{i}.m = inertial_true{i}.m * mass_scale(i);
        inertial_hat{i}.r_c = inertial_true{i}.r_c(:) * com_scale(i);
        I_hat = inertial_true{i}.I_c * inertia_scale(i);
        inertial_hat{i}.I_c = 0.5 * (I_hat + I_hat.');
    end

    info.mass_scale = mass_scale;
    info.com_scale = com_scale;
    info.inertia_scale = inertia_scale;
end
