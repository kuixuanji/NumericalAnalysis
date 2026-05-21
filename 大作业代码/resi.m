%% 残差函数
function[res, theta,P_e] = resi(psi, P_w, P_s, MDH_params, theta4_init, n_ref, k)
    a = MDH_params(:,1);
    alpha = MDH_params(:,2);
    d = MDH_params(:,3);
% 1.利用psi确定虚拟肘点坐标P_e
    D = norm(P_w - P_s);
    L_se = abs(d(3));%大臂长（不带偏置）
    L_ew = abs(sqrt(d(4)^2+d(5)^2));%小臂长（带偏置）
    % 上面如此定义的原因：由于腕点是确定的，L_EW不会因为关节5的旋转而改变,
    % 而L_se会因为关节3的旋转而改变，所以这里使用不带偏置的量，由于偏置引起
    % 的与实际肘点的偏移会在后面的迭代过程中被补偿

    %大臂向量与SW连线之间夹角的余弦值
    cos_a = (L_se^2+D^2-L_ew^2)/(2*L_se*D);
    % 边界安全检查，防止浮点数误差导致 acos 产生复数
    if cos_a > 1; cos_a = 1; end
    if cos_a < -1; cos_a = -1; end
    %虚拟肘点坐标P_e(使用几何推导，"手肘圆周"）
    P_e = P_s+L_se*cos_a*k+L_se*sqrt(1-cos_a^2)*(cos(psi)*cross(n_ref,k)+sin(psi)*n_ref);
% 2.反解前三个关节的旋转角度***
    theta(1) = -atan2(P_e(2),P_e(1))-asin(d(2)/sqrt(P_e(1)^2+P_e(2)^2));
    % 关节1下的肘点坐标(不等同于转移矩阵,只作平移与旋转theta1）
    xe_prime = P_e(1)*cos(theta(1))-P_e(2)*sin(theta(1));
    ye_prime = -P_e(1)*sin(theta(1))+P_e(2)*cos(theta(1));
    ze_prime = -P_e(3)-d(1);
    
    theta(2) = atan2(ye_prime, ze_prime);

    cos_theta3_component = ze_prime * cos(theta(2)) + ye_prime * sin(theta(2));
    theta(3) = atan2(d(2) - xe_prime, cos_theta3_component);
% 3.计算几何实际需要的 theta4***
    theta(4) = theta4_init;
    T_0_5 = eye(4);
    for i = 1:5
        if i == 5
        T_i = getMDH_Trans(a(i),alpha(i),d(i),0);
        else
        T_i = getMDH_Trans(a(i),alpha(i),d(i),theta(i));
        end
        T_0_5 = T_0_5*T_i;
    end
    P5_real = T_0_5(1:3, 4);

    res = norm(P_w - P5_real);
    
end