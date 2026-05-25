function Tau = RNEA(q,dq,ddq,MDH_Inertial)

%Recursive Newton-Euler Algorithm 递归牛顿欧拉算法

    % q, dq, ddq: 7x1 关节位置、速度、加速度向量
N = 7;
    % 预分配空间
    w = zeros(3, N);    % 角速度
    dw = zeros(3, N);   % 角加速度
    dv = zeros(3, N);   % 线加速度

    R_rel = zeros(3,3,N+1);% 存储连杆间的相对旋转矩阵
    P_rel = zeros(3,N+1);% 存储连杆间的相对平移向量
    
    % 初始条件
    w_prev = [0; 0; 0]; 
    dw_prev = [0; 0; 0];
    dv_prev = [0; 0; 9.81]; % 基座线加速度设为 9.81，自动计入重力补偿
    Z_axis = [0; 0; 1];     % 旋转轴均为 Z 轴
    
    MDH_params = [0,pi,-0.15643-0.12838;
                  0,pi/2,-0.005375-0.006375;
                  0,-pi/2,-0.21038-0.21038;
                  0,pi/2,-0.006375-0.006375;
                  0,-pi/2,-0.20843-0.10593;
                  0,pi/2,0;
                  0,-pi/2,-0.061525-0.10593;
                  0,pi,0];%%EE;
    %% Forward Pass
    for i = 1:size(MDH_params,1)
        a = MDH_params(i,1);
        alpha = MDH_params(i,2);
        d = MDH_params(i,3);
        if i == 8
            theta = 0;
        else
            theta = q(i);
        end
        tform = Get_MDH_Trans(a,alpha,d,theta);

        R_i_prev = tform(1:3,1:3);
        P_i_prev = tform(1:3, 4);

        R_rel(:,:,i) = R_i_prev;
        P_rel(:,i) = P_i_prev;

        % 获取转置矩阵：将 i-1 系下的向量投影到 i 系下
        R_prev_i = R_i_prev';

       % 只有前 7 个连杆需要计算运动学状态
        if i <= N
            % 1. 角速度
            w(:,i) = R_prev_i * w_prev + Z_axis * dq(i);

            % 2. 角加速度
            dw(:,i) = R_prev_i * dw_prev + cross(R_prev_i * w_prev, Z_axis * dq(i)) + Z_axis * ddq(i);

            % 3. 线加速度 (原点加速度)
            P_i = P_i_prev;
            dv(:,i) = R_prev_i * (dv_prev + cross(dw_prev, P_i) + cross(w_prev, cross(w_prev, P_i)));

            % 更新状态供下一步迭代
            w_prev = w(:,i);
            dw_prev = dw(:,i);
            dv_prev = dv(:,i);
        end
    end

    %% Backward Pass
    f = zeros(3, N+1); % 连杆间作用力
    n = zeros(3, N+1); % 连杆间力矩
    Tau = zeros(N, 1); % 输出关节力矩
    
    for i = N:-1:1
        m = MDH_Inertial{i}.m;
        rc = MDH_Inertial{i}.r_c;
        Ic = MDH_Inertial{i}.I_c;
        
        % 1. 计算当前连杆质心处的绝对线加速度 dv_c
        dv_c_i = dv(:, i) + cross(dw(:, i), rc) + cross(w(:, i), cross(w(:, i), rc));
        
        % 2. 计算维持当前运动所需的惯性力和力矩 (基于质心)
        Fi = m * dv_c_i; 
        Ni = Ic * dw(:,i) + cross(w(:,i), Ic * w(:,i));
        
        % 3. 获取下一级连杆传递回来的力和几何关系
            R_next = R_rel(:,:,i+1); % R_{i+1}^{i}
            p_next = P_rel(:,i+1);   % P_{i+1}^{i}
        
        % 4. 力的平衡方程 (牛顿第三定律，加上当前杆需要的力)
        f(:,i) = R_next * f(:,i+1) + Fi;
        
        % 5. 力矩的平衡方程 (包含上一级的力矩、力臂产生的力矩、自身的惯性力矩)
        n(:,i) = Ni + R_next * n(:,i+1) + cross(rc, Fi) + cross(p_next, R_next * f(:,i+1));
        
        % 6. 投影到当前关节的旋转轴 (mDH 的旋转轴始终是 Z 轴)
        Tau(i) = n(:,i)' * Z_axis;
    end
end
