function e_ori = Compute_QuatError(q_target, q_curr)
    wd = q_target(1); vd = q_target(2:4);
    wc = q_curr(1);   vc = q_curr(2:4);
    
    % delta_v 的哈密顿乘法展开
    delta_v = wd * (-vc) + wc * vd + cross(vd, -vc);
    
    e_ori = 2 * delta_v;
end