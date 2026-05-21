function theta = inverse_kinemics_hybrid(T_target,psi)

MDH_params = [0,pi,-0.1564-0.1284;
              0,pi/2,-0.0054-0.0064;
              0,-pi/2,-0.2104-0.2104;
              0,pi/2,-0.0064-0.0064;
              0,-pi/2,-0.2084-0.1059;
              0,pi/2,0;
              0,-pi/2,-0.0615-0.1059;
              0,pi,0];%%EE
a = MDH_params(:,1);
alpha = MDH_params(:,2);
d = MDH_params(:,3);
theta = zeros(7,1);
%% 求目标腕点坐标
T_7_EE = [1,0,0,0;
           0,-1,0,0;
           0,0,-1,0;
           0,0,0,1];
T_0_7 = T_target/T_7_EE;
P_w = T_0_7(1:3,4);
%% 猜测theta4，即迭代初值，（假设没有偏置的theta4）
L1 = abs(d(3));%大臂长
L2 = abs(d(5));%小臂长
P_s = [0;0;-d(1)];%肩点坐标
D = norm(P_w-P_s);%自运动轴线长度
cos_val = (L1^2 + L2^2 - D^2) / (2 * L1 * L2);

% 边界安全检查，防止浮点数误差导致 acos 产生复数
if cos_val > 1; cos_val = 1; end
if cos_val < -1; cos_val = -1; end
theta4_init = pi - acos(cos_val);
%% 建立固定参考平面
%自运动轴线向量k
k = (P_w-P_s)/norm(P_w-P_s);

%定义基座参考方向
Z0 = [0;0;1];
X0 = [1;0;0];
if abs(dot(k,Z0))>0.995 %点积绝对值接近1代表平行
    %此时发生奇异，k与Z0平行，改用X0轴构建参考平面
    n_ref = cross(k,X0)/norm(cross(k,X0));
    psi = psi-pi/2;%改变参考平面后需对psi进行修正
else %正常情况
    n_ref = cross(k,Z0)/norm(cross(k,Z0));
end
%% 一维数值迭代寻根
psi_now = psi;
max_iter = 15;
tol = 1e-7;
delta_psi = 1e-6;

tool = false;
theta_1_to_4 = zeros(4,1);

for iter = 1:max_iter
    % 1.计算当前残差
    [f_now,theta_now,P_e_now] = resi(psi_now, P_w, P_s, MDH_params, theta4_init, n_ref, k);
    
    if abs(f_now)<tol
        tool = true;
        theta_1_to_4 = theta_now;
        P_e = P_e_now;
        break;
    end

    % 2.计算微小扰动下的残差
    [f_pertur,~,~] = resi(psi_now+delta_psi,P_w,P_s,MDH_params,theta4_init,n_ref,k);

    % 3.计算一阶导数（斜率）
    df = (f_pertur - f_now) / delta_psi;

    % 避免分母过小导致数值发散（防奇异）
    if abs(df) < 1e-9
        df = sign(df) * 1e-9; 
    end
    
    % 4. 牛顿步进更新
    psi_now = psi_now - f_now / df;
end
% 迭代失败
if ~tool
        error('牛顿迭代法未收敛！可能目标位姿超出了工作空间，或者初始猜测值 psi 离解太远。');
end

% 迭代成功
theta(1:4) = theta_1_to_4;

%% 使用迭代结果求后三个关节角度
% 1.计算前四个关节的正向运动学矩阵
T_0_4 = eye(4);
for i = 1:4
    T_i = getMDH_Trans(a(i),alpha(i),d(i),theta(i));
    T_0_4 = T_0_4*T_i;
end
% 2.剥离出后三轴的正向运动学矩阵
T_4_7 = T_0_4\T_0_7;
% 3.提取左上角的 3x3 纯旋转矩阵 R_7_4
R_4_7 = T_4_7(1:3, 1:3);
% 4.解析反解后三个关节角度
r11 = R_4_7(1,1); r12 = R_4_7(1,2); r13 = R_4_7(1,3);
r21 = R_4_7(2,1); r22 = R_4_7(2,2); r23 = R_4_7(2,3);
r31 = R_4_7(3,1); r32 = R_4_7(3,2); r33 = R_4_7(3,3);
    % 奇异性阈值检查 (sin(theta6) 接近 0)
    sin_theta6 = sqrt(r31^2 + r32^2);
    
    if sin_theta6 < 1e-4
        % 【奇异状态处理】：J5 和 J7 共线
        theta(5) = 0; % 强行指定J5为0
        theta(6) = atan2(sin_theta6, r33); % 此时 theta6 接近 0 或 pi
        theta(7) = atan2(r21, r11);        % 全部旋转由 J7 承担
    else
        % 【正常状态处理】：根据矩阵元素对齐解析求解
        theta(6) = atan2(sin_theta6, r33);
        theta(5) = atan2(-r23, -r13);
        theta(7) = atan2(-r32, r31);
    end


end