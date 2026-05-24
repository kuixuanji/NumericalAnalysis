function tform = Forward_Kinematics(config)
    q = [config;0];
    MDH_params = [0,pi,-0.15643-0.12838;
              0,pi/2,-0.005375-0.006375;
              0,-pi/2,-0.21038-0.21038;
              0,pi/2,-0.006375-0.006375;
              0,-pi/2,-0.20843-0.10593;
              0,pi/2,0;
              0,-pi/2,-0.061525-0.10593;
              0,pi,0];%%EE
    tform = eye(4);
    for i = 1:size(MDH_params,1)
        a = MDH_params(i,1);
        alpha = MDH_params(i,2);
        d = MDH_params(i,3);
        theta = q(i);
        T_i = Get_MDH_Trans(a,alpha,d,theta);
        tform = tform*T_i;
    end
end