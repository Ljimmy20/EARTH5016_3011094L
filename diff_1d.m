% ***** 1D ADVECTION DIFFUSION MODEL OF HEAT TRANSPORT  *******************

% ***** Initialize Model Setup

% Set domain width and grid size
W = 16e3;    % Domain width [m]
Nx = 200;    % Number of grid cells in x-direction
h = W / Nx;  % Grid spacing based on domain width

% Create x-coordinate vectors
x_cells = h/2:h:W-h/2;  % Cell center positions in x-direction
x_faces = 0:h:W;        % Cell face positions in x-direction

% Set time step size
CFL = 0.25;              % Reduced CFL for more stability
u0 = 1;                 % Assumed constant velocity [m/s] (not used here)
k0 = 1;                 % Thermal diffusivity [m^2/s] (can be adjusted)
dt = CFL * min((h/2)/u0, (h/2)^2/k0);  % Time step based on CFL condition

% Set up index array for boundary conditions
BC = 'periodic';  % 'periodic' or 'insulating' based on the problem setup
switch BC
    case 'periodic'
        ind3 = [Nx, 1:Nx, 1];  % 3-point periodic boundary condition
        ind5 = [Nx-1, Nx, 1:Nx, 1, 2];  % 5-point periodic boundary condition
        ind7 = [Nx-2, Nx-1, Nx, 1:Nx, 1, 2, 3];  % 7-point periodic boundary condition
        
    case 'insulating'
        ind3 = [1, 1:Nx, Nx];  % 3-point insulating boundary condition
        ind5 = [1, 1, 1:Nx, Nx, Nx];  % 5-point insulating boundary condition
        ind7 = [1, 1, 1, 1:Nx, Nx, Nx, Nx];  % 7-point insulating boundary condition
end

% Set initial condition for temperature at cell centres
T0 = 20;            % Initial temperature [C]
dT = 10;            % Temperature change [C]
wT = 100;           % Width for initial Gaussian distribution [m]

T = T0 + dT .* exp(-(x_cells - W/2).^2 / (4 * wT^2));  % Gaussian initial temperature

% Store initial temperature and set initial temperature array
Tin = T;  % Store initial condition
Ta = T;   % Initialise analytical solution (for comparison later)

% *****  Solve Model Equations

% Time integration using Forward Euler (FE1)
t = 0;  % initial time [s]
tend = 16000;  % end time [s] (modify based on your needs)
output_interval = 10;  % Output interval in seconds (modify as needed)
k = 0;  % initial time step count

while t <= tend
    % Increment time and step count
    t = t + dt;
    k = k + 1;
    
    % Print current time for debugging
    disp(['Current time: ', num2str(t), ' s']);
    
    % Compute rate of change of temperature using diffusion only (no advection)
    dTdt = diffusion(T, k0, h, ind3);
    
    % Update temperature based on time step
    T = T + dTdt * dt;
    
    % Compute analytical solution at time t (for comparison)
    wTt = sqrt(wT^2 + 1 * k0 * t);
    Ta = T0 + dT .* wT ./ wTt .* exp(-(x_cells - W/2 - u0 * t).^2 / (4 * wTt^2)) + dT .* wT ./ wTt .* exp(-(x_cells + W/2 - u0 * t).^2 / (4 * wTt^2));
    
    % Calculate maximum temperature change
    temp_change = max(abs(T - Tin));  % Maximum temperature change
    
    % If the temperature change is small, stop the simulation
    if temp_change < 1e-6
        disp('Temperature change is too small, stopping the simulation.');
        break;
    end
    
    % Plot model progress
    if mod(k, output_interval) == 0
        makefig(x_cells, T, Tin, Ta, t);  % Example of plotting function
    end
end

% ***** Numerical error norm calculation
Err = norm(T - Ta, 2) / norm(Ta, 2);
disp(' ');
disp(['Numerical error = ', num2str(Err)]);
disp(' ');

% ***** Utility Functions  ************************************************

% Function to create output figure
function makefig(x, T, Tin, Ta, t)
    plot(x, Tin, 'k:', x, T, 'r-', x, Ta, 'k--', 'LineWidth', 1.5);
    axis tight; box on;
    xlabel('x [m]', 'FontSize', 15);
    ylabel('T [C]', 'FontSize', 15);
    title(['Temperature; time = ', num2str(t), ' s'], 'FontSize', 18);
    drawnow;
end

% Function to calculate diffusion rate
function dTdt = diffusion(f, k0, dx, ind)
    % Calculate heat flux by diffusion
    q = -k0 * diff(f(ind)) / dx;
    % Calculate flux balance for rate of change
    dTdt = -diff(q) / dx;
end