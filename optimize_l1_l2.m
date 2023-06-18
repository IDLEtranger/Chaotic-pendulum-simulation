clear; clc;
L1 = 0.15;
L2 = 0.1;
% 遍历法优化C、D点位置
t_max = 0;
l1_op = 0;
l2_op = 0;
re = 0;
% 优化 C 和 D 点位置
for l1_now = 0:L1/100:L1
    for l2_now = 0:L2/100:L2
        [t, r] = calculate(l1_now, l2_now);
        t_end = when_stop(l1_now, l2_now, r, t);
        if t_end >= t_max
        t_max = t_end;
        l1_op = l1_now;
        l2_op = l2_now;
        end
        re = re + 1;
    end
end

% 输出 t_max 的值
if t_max ~= 0
    fprintf("最长摆动时间：%f \n 最优l1、l2长度为：%f %f \n", t_max, l1_op, l2_op);
else
    fprintf("计算结果出错")
end