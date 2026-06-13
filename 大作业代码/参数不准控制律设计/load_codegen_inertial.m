function inertial = load_codegen_inertial(mat_file)
% 将代码生成用惯性参数文件读取为元胞数组格式。

    data = load(mat_file, 'm_all', 'rc_all', 'Ic_all');
    inertial = cell(1, 7);

    for i = 1:7
        inertial{i}.m = data.m_all(i);
        inertial{i}.r_c = data.rc_all(:, i);
        inertial{i}.I_c = data.Ic_all(:, :, i);
    end
end
