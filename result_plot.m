clear; clc;
% 常数定义
L1 = 0.15;
L2 = 0.1;
% l1 = 0.05;
% l2 = 0.025;
l1 = 0.040;
l2 = 0.089;

% D = 0.05; % DampingFactor Nxm/(rad/s)
D = 0.05;

% 计算混沌摆参数
tspan = [0 5];
r0 = [-pi/2 pi/2 0 0];
[t,r] = calculate(l1, l2, D, tspan, r0);

% 创建等间隔的时间点
N = ceil(length(t)/100)*100;
t_equ = linspace(tspan(1), tspan(2), N); % N为想要的时间点的数量

% 对 ODE 求解器的解进行插值
r_equ = interp1(t, r, t_equ);

% data
theta1 = r_equ(:,1);
theta2 = r_equ(:,2);
theta1_init = -pi/2;
theta2_init = pi;
point_C = [0, 0.5];

% 计算三角形的其他两个顶点坐标,用于后续绘制基点
edge_length = 0.05;
point_C1 = point_C + [0.5*edge_length*cos(pi/6), -edge_length*sin(pi/6)];
point_C2 = point_C + [0.5*edge_length*cos(5*pi/6), -edge_length*sin(5*pi/6)];

% equations of motion
point_D = [point_C(1) + l1*sin(theta1),...
           point_C(2) - l1*cos(theta1)];

point_E = [point_C(1) - (L1-l1)*sin(theta1),...
           point_C(2) + (L1-l1)*cos(theta1)];

point_F = [point_D(:,1) + (L2-l2)*sin(theta2),...
           point_D(:,2) - (L2-l2)*cos(theta2)];

point_G = [point_D(:,1) - l2*sin(theta2),...
           point_D(:,2) + l2*cos(theta2)];

% 绘制角度随时间变化图像
figure('Name','Solution of Lagrange function with ODE45');hold on;
plot(t_equ, r_equ(:,1), '-o', t_equ, r_equ(:,2), '-o', 'MarkerSize',3)
title('Solution of Lagrange function with ODE45');
xlabel('Time t');
ylabel('Solution q');
legend('q_1','q_2');

% 绘制E、F点位置随时间变化图像
figure('Name','x_Position of points E and F changing over time');hold on;
plot(t_equ,point_E(:,1),'-o',t_equ,point_F(:,1),'-o', 'MarkerSize',3)
title('x Position of points E and F');
xlabel('Time t');
ylabel('Position x');
legend('E_x','F_x');

figure('Name','y_Position of points E and F changing over time');hold on;
plot(t_equ,point_E(:,2),'-o',t_equ,point_F(:,2),'-o', 'MarkerSize',3)
title('y Position of points E and F');
xlabel('Time t');
ylabel('Position x');
legend('E_y','F_y');
%%%%%%%%%%%%%%%%%%%%%%%
% 计算 E 点和 F 点的速度和加速度
dt = t_equ(2) - t_equ(1);  % 时间间隔

% 计算 E 点的速度和加速度
velocity_E_x = gradient(point_E(:,1)) / dt;
velocity_E_y = gradient(point_E(:,2)) / dt;
acceleration_E_x = gradient(velocity_E_x) / dt;
acceleration_E_y = gradient(velocity_E_y) / dt;

% 计算 F 点的速度和加速度
velocity_F_x = gradient(point_F(:,1)) / dt;
velocity_F_y = gradient(point_F(:,2)) / dt;
acceleration_F_x = gradient(velocity_F_x) / dt;
acceleration_F_y = gradient(velocity_F_y) / dt;

% 绘制速度图像
figure('Name','x Velocity of points E and F'); hold on;
plot(t_equ, velocity_E_x, '-o', t_equ, velocity_F_x, '-o', 'MarkerSize',3)
title('x Velocity of points E and F');
xlabel('Time t');
ylabel('Velocity_x');
legend('v_Ex','v_Fx');

figure('Name','y Velocity of points E and F'); hold on;
plot(t_equ, velocity_E_y, '-o', t_equ, velocity_F_y, '-o', 'MarkerSize',3)
title('y Velocity of points E and F');
xlabel('Time t');
ylabel('Velocity_y');
legend('v_Ey','v_Fy');

% 绘制加速度图像
figure('Name','x Acceleration of points E and F'); hold on;
plot(t_equ, acceleration_E_x, '-o', t_equ, acceleration_F_x, '-o', 'MarkerSize',3)
title('x Acceleration of points E and F');
xlabel('Time t');
ylabel('Acceleration_x');
legend('a_Ex','a_Fx');
figure('Name','y Acceleration of points E and F'); hold on;
plot(t_equ, acceleration_E_y, '-o', t_equ, acceleration_F_y, '-o', 'MarkerSize',3)
title('y Acceleration of points E and F');
xlabel('Time t');
ylabel('Acceleration_y');
legend('a_Ey','a_Fy');

% 提取所有点的x和y坐标
xD = point_D(:,1);
yD = point_D(:,2);

xE = point_E(:,1);
yE = point_E(:,2);

xF = point_F(:,1);
yF = point_F(:,2);

xG = point_G(:,1);
yG = point_G(:,2);

% 计算所有点的x和y的范围
xMin = min([xD; xE; xF; xG]) - 0.1;
xMax = max([xD; xE; xF; xG]) + 0.1;
yMin = min([yD; yE; yF; yG]) - 0.1;
yMax = max([yD; yE; yF; yG]) + 0.1;
width = max(xMax - xMin, yMax - yMin);

% 创建一个VideoWriter对象
v = VideoWriter(sprintf('Simulation_l1(%.3f)_l2(%.3f)_D(%.3f).avi', l1, l2, D));
v.Quality = 90;  % 设置视频质量为 90
% 设置帧率为1/dt帧每秒
v.FrameRate = length(t_equ)/tspan(2);
% 打开视频文件
open(v);
figure('Name','Linkage Simulation'); hold on;
for k = 1:length(t_equ)
    % 清空当前的图像，准备画新的一帧
    cla;
    % 绘制基点
    % 绘制等边三角形
    vertices = [point_C1; point_C2; point_C];
    faces = [1, 2, 3, 1];
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', [0.3010 0.7450 0.9330], 'EdgeColor', [0 0.4470 0.7410]);

    % 在右上角显示当前的时间步长
    % 在右上角显示当前的时间步长，保留三位小数，并加上前缀't = '
    text(xMax, yMax, sprintf('t = %.3f', t_equ(k)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    
    % 设置坐标轴的范围，使C点在视频中位置保持静止，xy轴比例
    xlim([point_C(1) - width/2, point_C(1) + width/2]);
    ylim([point_C(2) - width/2, point_C(2) + width/2]);
    pbaspect([1 1 1]);
    
    % 画出当前位置
    plot(xD(k), yD(k), '.', 'MarkerSize',20);
    plot(xE(k), yE(k), '.', 'MarkerSize',20);
    plot(xF(k), yF(k), '.', 'MarkerSize',20);
    plot(xG(k), yG(k), '.', 'MarkerSize',20);
    plot(point_C(1), point_C(2), '.', 'MarkerSize',20);

    % 连接点E和D，以及点G和F
    line([xE(k) xD(k)], [yE(k) yD(k)], 'LineWidth',1 , 'Color', 'k');
    line([xG(k) xF(k)], [yG(k) yF(k)], 'LineWidth',1 , 'Color', 'k');

    % 获取当前帧
    frame = getframe;
    
    % 将帧写入视频
    writeVideo(v, frame);
end
% 关闭视频文件
close(v);

filePath = mfilename('fullpath');
[path, ~, ~] = fileparts(filePath);
disp("写入完成！视频文件存放于：")
disp(path);