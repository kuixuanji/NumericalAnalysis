clear; clc;

load('MDH_Inertial.mat', 'MDH_Inertial');

m_all  = zeros(7, 1);
rc_all = zeros(3, 7);
Ic_all = zeros(3, 3, 7);

for i = 1:7
    m_all(i) = MDH_Inertial{i}.m;
    rc_all(:, i) = MDH_Inertial{i}.r_c(:);
    Ic_all(:, :, i) = MDH_Inertial{i}.I_c;
end

save('MDH_Inertial_Codegen.mat', 'm_all', 'rc_all', 'Ic_all');