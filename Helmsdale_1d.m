%***** 1D HEAT MODEL BASED ON HELMSDALE DATA *********************************

% clear workspace
clear all; close all; clc;

%***** PARAMETERS ***********************************************************
% Load data directly (replace these with values extracted from the Excel)
depth = [51, 150, 249, 351]; % Depth in meters
temperature = [6.5, 11.4, 15.7, 21.1]; % Temperature in Celsius
conductivity = [0.6239, 0.7429, 0.7853, 0.9103]; % Thermal conductivity [W/m/K]
density = [1837, 1861, 1880, 1894]; % Density [kg/m^3]
heat_capacity = 1000; % Heat capacity [J/kg/K] (assumed uniform)

% Domain setup
W = max(depth); % Total depth [m]
Nx = 200; % Number of grid points
dz = W / Nx; % Grid spacing

% Time parameters
yr = 3600 * 24 * 365; % Seconds in a year
tend = 1e5 * yr; % Simulation time [s]
CFL = 0.5; % CFL number for stability

% Thermal diffusivity (kT / (rho * Cp))
thermal_diffusivity = conductivity ./ (density * heat_capacity);

% Interpolate depth-dependent properties to grid
z = linspace(0, W, Nx);
kT = interp1(depth, conductivity, z, 'linear', 'extrap');
rho = interp1(depth, density, z, 'linear', 'extrap');
Cp = heat_capacity * ones(size(z)); % Assumed constant

% Initial temperature profile (linear gradient)
T = interp1(depth, temperature, z, 'linear', 'extrap');

% Boundary conditions
T_top = T(1); % Surface temperature
T_bottom = T(end); % Bottom temperature

% Time step calculation
dt = CFL * dz^2 / max(thermal_diffusivity);

%***** SIMULATION ***********************************************************

t = 0; % Initial time
while t < tend
    % Ensure T is a column vector
    T = T(:); % Convert to column vector

    % Compute heat flux by diffusion at cell faces
    T_extended = [T_top; T; T_bottom]; % Extend temperature with boundary conditions
    
    % Calculate q at internal cell faces, it will have size Nx-1
    q = -diff(kT .* T_extended) / dz; % Heat flux at cell faces

    % Initialize q_extended with correct size Nx
    q_extended = zeros(Nx, 1); % Initialize with zeros, length Nx
    
    % Assign the internal fluxes from q (q has size Nx-1)
    q_extended(2:end-1) = q; % Corrected to match size with q (Nx-1)

    % Apply boundary conditions for q
    % Here, we set the boundary flux to zero or some other fixed value as needed
    q_extended(1) = 0; % Boundary at the top
    q_extended(end) = 0; % Boundary at the bottom
    
    % Initialize dTdt array
    dTdt = zeros(Nx, 1); 

    % Compute dTdt for interior cells (excluding boundary)
    dTdt(2:end-1) = -(q_extended(2:end) - q_extended(1:end-1)) ./ (rho(2:end-1) .* Cp(2:end-1)); % Corrected for interior cells
    
    % Update temperature (forward Euler)
    T = T + dTdt * dt;

    % Increment time
    t = t + dt;
end

%***** PLOT RESULTS *********************************************************
figure;
plot(z, T, 'r-', 'LineWidth', 2);
xlabel('Depth [m]', 'FontSize', 12);
ylabel('Temperature [°C]', 'FontSize', 12);
title('1D Temperature Distribution', 'FontSize', 14);
grid on;
