function t_end = when_stop(l1, l2, r, t)
% 判断系统是否达到了设定的终止条件
    % 常数定义
    L1 = 0.15;
    L2 = 0.1;
    %l1 = 0.05;
    %l2 = 0.025;
    
    m1 = 0.6;
    m2 = 0.5;
    
    J1 = (1/12) * (m1*L1^2); % moment of inertia
    J2 = (1/12) * (m2*L2^2); 
    
    D = 0.05; % DampingFactor Nxm/(rad/s)
    g = 9.80; % Gravitational acceleration

    % 系统最小机械能
    q1_end = r(end, 1);
    q2_end = r(end, 2);
    % 定义要圆整的幅度（正负10度对应的弧度）
    tolerance = 10 * pi / 180;

    % 圆整 q_end 到最近的 π 整倍数
    rounded_q1_end = round(q1_end / pi) * pi;
    rounded_q2_end = round(q2_end / pi) * pi;
    
    % 检查圆整后的值是否超出指定的容差范围
    if (abs(q1_end - rounded_q1_end) < tolerance) && (abs(q2_end - rounded_q2_end) < tolerance)
        q1_end = rounded_q1_end;
        q2_end = rounded_q2_end;
    else
        fprintf('warning：tspan时间太短'); % 机构在超过tspan时仍未趋近稳态
    end

    Ek_min = 0;
    Ep_min = m1*g*(0.5 + (0.5*L1-l1)*cos(q1_end)) + m2*g*(0.5 - (l1*cos(q1_end)) - (0.5*L2 - l2)*cos(q2_end));
    Emin = Ep_min + Ek_min;
    
    for n = 1:size(r,1)
        % 系统机械能计算
        q1 = r(n, 1);
        q2 = r(n, 2);
        dq1 = r(n, 3);
        dq2 = r(n, 4);
        % 偏类速度
        ux1_1 = -(0.5*L1 - l1)*cos(q1);
        ux1_2 = 0;
        uy1_1 = -(0.5*L1 - l1)*sin(q1);
        uy1_2 = 0;
        utheta1_1 = 1;
        utheta1_2 = 0;
        
        ux2_1 = l1*cos(q1);
        ux2_2 = (0.5*L2 - l2)*cos(q2);
        uy2_1 = l1*sin(q1);
        uy2_2 = (0.5*L2 - l2)*sin(q2);
        utheta2_1 = 0;
        utheta2_2 = 1;
        
        % 等效转动惯量
        J11 = m1*((ux1_1)^2 + (uy1_1)^2) + J1*(utheta1_1)^2 + m2*((ux2_1)^2 + (uy2_1)^2) + J2*(utheta2_1)^2;
        J22 = m1*((ux1_2)^2 + (uy1_2)^2) + J1*(utheta1_2)^2 + m2*((ux2_2)^2 + (uy2_2)^2) + J2*(utheta2_2)^2;
        J12 = m1*(ux1_1*ux1_2 + uy1_1*uy1_2) + J1*utheta1_1*utheta1_2 + m2*(ux2_1*ux2_2 + uy2_1 * uy2_2) + J2*utheta2_1*utheta2_2;
    
        % 系统机械能
        Ek = (1/2)*J11*dq1^2 + (1/2)*J22*dq2^2 + J12*dq1*dq2;
        Ep = m1*g*(0.5 + (0.5*L1-l1)*cos(q1)) + m2*g*(0.5 - (l1*cos(q1)) - (0.5*L2 - l2)*cos(q2));
        E = Ep + Ek;
    
        if E <= 1.001*Emin
            t_end = t(n);
            break
        end
    end
    
end
