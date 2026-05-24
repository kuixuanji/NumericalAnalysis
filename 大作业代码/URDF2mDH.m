clc;clear;close all;    
%这个脚本得到：
%-mDH参数描述坐标系下的质心坐标
%-参考点为mDH坐标系原点且轴向与mDH坐标系的各轴平行的惯性张量


% 初始化 7 连杆的惯性参数 Cell 数组
% 各字段定义：
% m: 质量 (kg)
% r_c: 质心向量 (3x1, 相对于URDF各连杆自身坐标系)
% I_c: 惯性张量 (3x3, 相对于质心系)
URDF_Inertial = cell(7, 1);

% --- Link 1: shoulder_link ---
URDF_Inertial{1}.m = 1.3773;
URDF_Inertial{1}.r_c = [-2.3E-05; -0.010364; -0.07336];
URDF_Inertial{1}.I_c = [0.00457, 1E-06, 2E-06; 
                        1E-06, 0.004831, 0.000448; 
                        2E-06, 0.000448, 0.001409];

% --- Link 2: half_arm_1_link ---
URDF_Inertial{2}.m = 1.1636;
URDF_Inertial{2}.r_c = [-4.4E-05; -0.09958; -0.013278];
URDF_Inertial{2}.I_c = [0.011088, 5E-06, 0; 
                        5E-06, 0.001072, -0.000691; 
                        0, -0.000691, 0.011255];

% --- Link 3: half_arm_2_link ---
URDF_Inertial{3}.m = 1.1636;
URDF_Inertial{3}.r_c = [-4.4E-05; -0.006641; -0.117892];
URDF_Inertial{3}.I_c = [0.010932, 0, -7E-06; 
                        0, 0.011127, 0.000606; 
                        -7E-06, 0.000606, 0.001043];

% --- Link 4: forearm_link ---
URDF_Inertial{4}.m = 0.9302;
URDF_Inertial{4}.r_c = [-1.8E-05; -0.075478; -0.015006];
URDF_Inertial{4}.I_c = [0.008147, -1E-06, 0; 
                        -1E-06, 0.000631, -0.0005; 
                        0, -0.0005, 0.008316];

% --- Link 5: spherical_wrist_1_link ---
URDF_Inertial{5}.m = 0.6781;
URDF_Inertial{5}.r_c = [1E-06; -0.009432; -0.063883];
URDF_Inertial{5}.I_c = [0.001596, 0, 0; 
                        0, 0.001607, 0.000256; 
                        0, 0.000256, 0.000399];

% --- Link 6: spherical_wrist_2_link ---
URDF_Inertial{6}.m = 0.6781;
URDF_Inertial{6}.r_c = [1E-06; -0.045483; -0.00965];
URDF_Inertial{6}.I_c = [0.001641, 0, 0; 
                        0, 0.00041, -0.000278; 
                        0, -0.000278, 0.001641];

% --- Link 7: bracelet_link ---
URDF_Inertial{7}.m = 0.5006;
URDF_Inertial{7}.r_c = [-0.000281; -0.011402; -0.029798];
URDF_Inertial{7}.I_c = [0.000587, 3E-06, 3E-06; 
                        3E-06, 0.000369, 0.000118; 
                        3E-06, 0.000118, 0.000609];
MDH_Inertial = cell(7, 1);
for i = 1:7
    MDH_Inertial{i}.m = URDF_Inertial{i}.m;
    %获取转移矩阵(Link i 的 URDF系 到 mDH系)
    T_bridge = get_mDH2URDF(i);
    R_bridge = T_bridge(1:3,1:3);
    p_bridge = T_bridge(1:3,4);
    % 1.质心转换
    MDH_Inertial{i}.r_c = R_bridge*URDF_Inertial{i}.r_c+p_bridge;

    % 2.惯性张量转换 (旋转 + 平行轴定理)
    Ic_urdf = URDF_Inertial{i}.I_c;
    
    % 先旋转
    Ic_rotated = R_bridge * Ic_urdf * R_bridge';
    
    % 再使用平行轴定理 (Steiner's Theorem)
    % 公式: I_new = I_rot + m * ( (r'*r)*E - r*r' )
    r = MDH_Inertial{i}.r_c;
    Steiner = MDH_Inertial{i}.m * ( (r'*r)*eye(3) - (r*r') );
    
    MDH_Inertial{i}.I_c = Ic_rotated;
end
save('MDH_Inertial.mat','MDH_Inertial');


function T = get_mDH2URDF(i)
    T_DH = get_mDH_transform(i);
    T_URDF = get_urdf_transform(i);
    T = T_DH\T_URDF;                                                                               
end

function T = get_mDH_transform(i)
        MDH_params = [0,pi,-0.15643-0.12838;
              0,pi/2,-0.005375-0.006375;
              0,-pi/2,-0.21038-0.21038;
              0,pi/2,-0.006375-0.006375;
              0,-pi/2,-0.20843-0.10593;
              0,pi/2,0;
              0,-pi/2,-0.061525-0.10593;
              0,pi,0];
        T = eye(4);
        for j = 1:i
        a = MDH_params(j,1);
        alpha = MDH_params(j,2);
        d = MDH_params(j,3);
        theta = 0;
        T_j = Get_MDH_Trans(a,alpha,d,theta);
        T = T*T_j;
        end
end

function T = get_urdf_transform(i)
    T = eye(4);
    for j = 1:i
        T_now = get_joint_transform(j);
        T = T*T_now;
    end
end

function T = get_joint_transform(index)
    % 1. 定义从 URDF 中提取的零位常数偏移矩阵
    % 数据已在上文中结构化展开，这里写成紧凑的逻辑
    T_offset = eye(4);
    switch index
        case 1
            T_offset = [1,0,0,0; 0,-1,0,0; 0,0,-1,0.15643; 0,0,0,1];
        case 2
            T_offset = [1,0,0,0; 0,0,-1,0.005375; 0,1,0,-0.12838; 0,0,0,1];
        case 3
            T_offset = [1,0,0,0; 0,0,1,-0.21038; 0,-1,0,-0.006375; 0,0,0,1];
        case 4
            T_offset = [1,0,0,0; 0,0,-1,0.006375; 0,1,0,-0.21038; 0,0,0,1];
        case 5
            T_offset = [1,0,0,0; 0,0,1,-0.20843; 0,-1,0,-0.006375; 0,0,0,1];
        case 6
            T_offset = [1,0,0,0; 0,0,-1,0.00017505; 0,1,0,-0.10593; 0,0,0,1];
        case 7
            T_offset = [1,0,0,0; 0,0,1,-0.10593; 0,-1,0,-0.00017505; 0,0,0,1];
    end
       
    % 3. 复合得到最终的实时齐次变换矩阵
    T = T_offset;
end
