function drdt = chaos(t, r, ddq1_str, ddq2_str)
% r(1) = q1 r(2) = q2 r(3) = dq1 r(4) = dq2
    drdt = zeros(4, 1);
    drdt(1) = r(3); % r1' = r3
    drdt(2) = r(4); % r2' = r4
    drdt(3) = eval(ddq1_str); % r3' = ddq1
    drdt(4) = eval(ddq2_str); % r4' = ddq2
end