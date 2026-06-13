function y = sat_vector(x)
% 对向量逐元素限幅到 [-1, 1]。
% 饱和函数
    y = min(max(x, -1), 1);
end
