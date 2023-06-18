function [t, r] = calculate(l1, l2, D, tspan, r0)
    % 计算混沌摆状态
    if nargin < 1
        l1 = 0.05;
    end
    if nargin < 2
        l2 = 0.025;
    end
    if nargin < 3
        D = 0.05; % DampingFactor Nxm/(rad/s)
    end
    if nargin < 4
        tspan = [0 1000];
    end
    if nargin < 5
        r0 = [-pi/2 pi/2 0 0];
    end
    % 系统中常量
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    L1 = 0.15;
    L2 = 0.1;

    m1 = 0.6;
    m2 = 0.5;   

    J1 = (1/12) * (m1*L1^2); % moment of inertia
    J2 = (1/12) * (m2*L2^2); 
    
    %D = 0.05; % DampingFactor Nxm/(rad/s)
    g = 9.80; % Gravitational acceleration
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % 广义坐标及其导数
    syms q1(t) q2(t) 
    dq1 = diff(q1, t);
    dq2 = diff(q2, t);
    ddq1 = diff(dq1, t);
    ddq2 = diff(dq2, t);
    
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
    
    %syms F1y F2y MC MD % 外力、外力矩
    F1y = -m1*g;
    F2y = -m2*g;
    MC = -D*dq1;
    MD = -D*dq2;
    M1 = MC - MD;
    M2 = MD;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    Ek = (1/2)*J11*dq1^2 + (1/2)*J22*dq2^2 + J12*dq1*dq2;
    Ep = m1*g*(0.5 + (0.5*L1-l1)*cos(q1)) + m2*g*(0.5 - (l1*cos(q1)) - (0.5*L2 - l2)*cos(q2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q1 = F1y*uy1_1 + F2y*uy2_1 + M1*utheta1_1 + M2*utheta2_1;
    Q2 = F1y*uy1_2 + F2y*uy2_2 + M1*utheta1_2 + M2*utheta2_2;
    
    %%%%%%%%%%%%% 求出将要带入lagrange方程的项 %%%%%%%%%%%%
    roundEk_round_dq1 = diff(Ek, dq1);
    d_roundEk_round_dq1 = diff(roundEk_round_dq1, t);
    roundEk_roundq1 = diff(Ek, q1);
    roundEp_roundq1 = diff(Ep, q1);
    
    roundEk_round_dq2 = diff(Ek, dq2);
    d_roundEk_round_dq2 = diff(roundEk_round_dq2, t);
    roundEk_roundq2 = diff(Ek, q2);
    roundEp_roundq2 = diff(Ep, q2);
    
    %%%%%%%%%%%%%%%%% 列写lagrange方程 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eqn1 = d_roundEk_round_dq1 - roundEk_roundq1 + roundEp_roundq1 == Q1;
    eqn2 = d_roundEk_round_dq2 - roundEk_roundq2 + roundEp_roundq2 == Q2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    syms r_tmp [1 6]
    eqn1_s = subs(eqn1, [q1, q2, dq1, dq2, ddq1, ddq2], r_tmp);
    eqn2_s = subs(eqn2, [q1, q2, dq1, dq2, ddq1, ddq2], r_tmp);
    [ddq1_sol, ddq2_sol] = solve([eqn1_s, eqn2_s], [r_tmp(5), r_tmp(6)]);    

    ddq1_str = char(ddq1_sol);
    ddq2_str = char(ddq2_sol);
    % 对字符串进行替换
    ddq1_str = strrep(ddq1_str, 'r_tmp1', 'r(1)');
    ddq1_str = strrep(ddq1_str, 'r_tmp2', 'r(2)');
    ddq1_str = strrep(ddq1_str, 'r_tmp3', 'r(3)');
    ddq1_str = strrep(ddq1_str, 'r_tmp4', 'r(4)');
    
    ddq2_str = strrep(ddq2_str, 'r_tmp1', 'r(1)');
    ddq2_str = strrep(ddq2_str, 'r_tmp2', 'r(2)');
    ddq2_str = strrep(ddq2_str, 'r_tmp3', 'r(3)');
    ddq2_str = strrep(ddq2_str, 'r_tmp4', 'r(4)');

    % 求解
    % tspan = [0 5];
    % r0 = [-pi/2 pi/2 0 0];
    % 匿名函数传递参数给 chaos 函数
    odefun = @(t, r) chaos(t, r, ddq1_str, ddq2_str);
    % 求解微分方程
    [t, r] = ode23(odefun, tspan, r0);
end