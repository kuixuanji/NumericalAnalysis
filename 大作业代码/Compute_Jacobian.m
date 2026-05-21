function J = Compute_Jacobian(q,MDH_params)
    % q: 当前 7个关节角 [q1, q2, ..., q7]
    num_joints = 7;
    J = zeros(6,num_joints);

    T_cell = cell(1,num_joints+1);
    T_current = eye(4);
    for i = 1:num_joints
        a = MDH_params(i,1);
        alpha = MDH_params(i,2);
        d = MDH_params(i,3);
        theta = q(i);

        T_mid = Get_MDH_Trans(a,alpha,d,theta);
        T_current = T_current*T_mid;
        T_cell{i} = T_current;
    end
    % 加上最后一行到末端手爪（Tool）的固定变换
    if size(MDH_params, 1) > num_joints
        a = MDH_params(8, 1); alpha = MDH_params(8, 2); d = MDH_params(8, 3);
        T_tool = Get_MDH_Trans(a,alpha,d,0);
        T_current = T_current * T_tool;
    end
    
    p_n = T_current(1:3,4);%末端手爪在基系下的位置坐标

    %构造雅可比矩阵
    for i = 1:num_joints
        T_i = T_cell{i};

        z_i = T_i(1:3, 3); % 提取Z轴方向
        p_i = T_i(1:3, 4); % 提取原点位置

        % 叉乘计算线速度贡献，直接赋值角速度贡献
        J_v = cross(z_i, (p_n - p_i));
        J_w = z_i;

        J(:,i) = [J_v;J_w];
    end
end
