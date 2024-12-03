%***** RUN 1D ADVECTION DIFFUSION MODEL WITH BOREHOLE LOG DATA ***********************************

% clear workspace
clear all; close all; %clc;

% 读取孔隙数据
data = readtable('Borehole.csv'); % 假设数据保存在名为 'borehole_data.csv' 的文件中

% set model parameters based on borehole data
D = data.D;      % 深度 [m]
Type = data.Type; % 材料类型
A = data.A;       % 年龄 [Myr]
T = data.T;       % 温度 [C]
kT = data.kT;     % 热导率 [W/M/K]
rho = data.rho;   % 密度 [kg/m3]

% Define grid size and spacing
N = length(D);    % 网格数目
dx = D(2) - D(1); % 假设深度间隔一致

% Create coordinate vectors for grid cells and faces
xc = D; % 使用深度数据作为坐标
xf = [0; D]; % 增加边界点

% Time step size
CFL = 1/16;         % Time step limiter
dt = CFL * min((dx/2) / 1e-6, (dx/2)^2 / 1e-6); % Time step [s]

% Set initial condition for temperature at cell centers
T0 = T;             % 使用温度数据作为初始条件
Tin = T0;           % 存储初始条件
Ta = T0;            % 初始的解析解

% Set up boundary conditions (using periodic example)
BC = 'periodic';    % Boundary condition option flag ('insulating', 'periodic')
switch BC
    case 'periodic'
        ind3 = [N, 1:N, 1];
        ind5 = [N-1, N, 1:N, 1, 2];
        ind7 = [N-2, N-1, N, 1:N, 1, 2, 3];
    case 'insulating'
        ind3 = [1, 1:N, N];
        ind5 = [1, 1, 1:N, N, N];
        ind7 = [1, 1, 1, 1:N, N, N, N];
end

% Initialize output figure with initial condition
figure(1); clf
makefig(xc, T0, Tin, Ta, 0);

%***** RUN MODEL
tend = D(end)/1e-6;   % Ending time [s] based on grid size and advection speed
t = 0;                % Initial time [s]
k = 0;                % Initial time step count

%***** Solve Model Equations
while t <= tend
    % Increment time and step count
    t = t + dt;
    k = k + 1;

    % Apply the Forward Euler time integration scheme
    dTdt = diffusion(T, kT, dx, ind3) + advection(T, 1e-6, dx, ind7, 'WENO5', kT, rho);

    T = T + dTdt * dt;  % Update temperature field

    % Get analytical solution at time t (this can be customized based on your specific solution)
    Ta = T0; % For simplicity, assume analytical solution is equal to initial condition

    % Plot model progress every 'nop' steps
    if mod(k, 100) == 0
        makefig(xc, T, Tin, Ta, t / (3600 * 24 * 365)); % Convert time to years for plot
    end
end

% Calculate numerical error norm
Err = norm(T - Ta, 2) / norm(Ta, 2);
disp(' ');
disp(['Numerical error = ', num2str(Err)]);
disp(' ');

%***** Utility Functions

% Function to make output figure
function makefig(x, T, Tin, Ta, t)
    plot(x, Tin, 'k:', x, T, 'r-', x, Ta, 'k--', 'LineWidth', 1.5); axis tight; box on;
    xlabel('Depth [m]', 'FontSize', 15);
    ylabel('Temperature [°C]', 'FontSize', 15);
    title(['Temperature; time = ', num2str(t), ' years'], 'FontSize', 18);
    drawnow;
end

% Function to calculate diffusion rate
function dTdt = diffusion(f, kT, dx, ind)
    % Calculate heat flux by diffusion
    q = - kT .* diff(f(ind)) / dx;
    % Calculate flux balance for rate of change
    dTdt = - diff(q) / dx;
end

% Function to calculate advection rate
function dTdt = advection(f, u0, dx, ind, ADVN, kT, rho)
    % Split velocities into positive and negative
    upos = 0.5 * (u0 + abs(u0));    % Positive velocity
    uneg = 0.5 * (u0 - abs(u0));    % Negative velocity

    % Get values on stencil nodes
    fmmm = f(ind(1:end-6));
    fmm  = f(ind(2:end-5));
    fm   = f(ind(3:end-4)); 
    fc   = f(ind(4:end-3)); 
    fp   = f(ind(5:end-2)); 
    fpp  = f(ind(6:end-1));
    fppp = f(ind(7:end-0));

    % Calculate heat flux by advection
    switch ADVN
        case 'UPW1'
            fppos = fc;     fpneg = fp;
            fmpos = fm;     fmneg = fc;

        case 'WENO5'
            fppos = weno5poly(fmm ,  fm, fc, fp, fpp); 
            fpneg = weno5poly(fppp, fpp, fp, fc, fm );
            fmpos = weno5poly(fmmm, fmm, fm, fc, fp ); 
            fmneg = weno5poly(fpp, fp, fc, fm, fmm);
    end

    % Calculate flux balance for rate of change
    div_qpos = upos .* (fppos - fmpos) / dx;
    div_qneg = uneg .* (fpneg - fmneg) / dx;
    div_q    = div_qpos + div_qneg;

    dTdt = - div_q;
end

% WENO5 interpolation function
function [fface] = weno5poly(fmm, fm, fc, fp, fpp)
    % 5th order WENO polynomials
    p1 = (2*fmm - 7*fm + 11*fc ) / 6;
    p2 = ( -fm  + 5*fc + 2*fp ) / 6;
    p3 = (2*fc  + 5*fp - fpp) / 6;

    % Smoothness measure
    b1 = 13/12 * (fmm - 2*fm + fc).^2 + 1/4 * (fmm - 4*fm + 3*fc).^2;
    b2 = 13/12 * (fm  - 2*fc + fp).^2 + 1/4 * (fm  - fp).^2;
    b3 = 13/12 * (fc  - 2*fp + fpp).^2 + 1/4 * (3*fc - 4*fp + fpp).^2;

    % Weights
    g   = [1/10, 6/10, 3/10];
    wp1 = g(1) ./ (b1.^2 + eps);
    wp2 = g(2) ./ (b2.^2 + eps);
    wp3 = g(3) ./ (b3.^2 + eps);

    % Get flux (normalize weights at the same time)
    fface = (wp1 .* p1 + wp2 .* p2 + wp3 .* p3) ./ (wp1 + wp2 + wp3);
end

