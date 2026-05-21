function tform = Forward_Kinematics(config)
    q = [config;0];
    MDH_params = [0,pi,-0.1564-0.1284;
                  0,pi/2,-0.0054-0.0064;
                  0,-pi/2,-0.2104-0.2104;
                  0,pi/2,-0.0064-0.0064;
                  0,-pi/2,-0.2084-0.1059;
                  0,pi/2,0;
                  0,-pi/2,-0.0615-0.1059;
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