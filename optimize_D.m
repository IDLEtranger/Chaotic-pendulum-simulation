clear;clc;
D_max = 1;
t_max = 0;
D_op = 1;
l1_op = 0.0405;
l2_op = 0.0890;
flag = 1;
D = zeros(100, 1);
t_D = zeros(100, 1);
% 遍历法优化D
% 优化 C 和 D 点位置
for D_now = D_max:-D_max/100:0
    [t, r] = calculate(l1_op, l2_op, D_now);
    t_end = when_stop(l1_op, l2_op, r, t);
    D(flag) = D_now;
    t_D(flag) = t_end;
    if t_end >= t_max
        t_max = t_end;
        D_op = D_now;
    end
    flag = flag + 1;
end

% 输出 t_max 的值
if t_max ~= 0
    fprintf("最长摆动时间：%f \n 最优D为：%f \n", t_max, D_op);
else
    fprintf("计算结果出错")
end