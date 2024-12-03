%*****  Initialize Model Setup  ****************************************

% Parameters
W = 1000;   % width of the domain [m]
D = 500;    % depth of the domain [m]
h = 10;     % grid spacing [m]
Ttop = 0;   % top temperature [°C]
Tbot = 100; % bottom temperature [°C]
kT = 2;     % thermal conductivity [W/mK]
alpha = kT / (1000); % Thermal diffusivity (assuming density and specific heat)
dt_yr = 1;  % time step [years]
tend = 100; % total time to simulate [years]

% create x-coordinate vectors
xc = h/2:h:W-h/2;      % x-coordinate vector for cell centre positions [m]
zc = h/2:h:D-h/2;      % z-coordinate vector for cell centre positions [m]
[Xc,Zc] = meshgrid(xc,zc);  % create 2D coordinate arrays

% set initial condition for temperature at cell centres (linear gradient)
T = Ttop + (Tbot - Ttop) * Zc / D;  

% plot initial condition
figure;
h_img = imagesc(xc, zc, T); 
axis equal tight;
colorbar;
caxis([Ttop Tbot]);
title(['Temperature Distribution at t = ', num2str(0, '%.1f'), ' years']);
xlabel('Distance (m)');
ylabel('Depth (m)');
drawnow;

% Add colorbar once and use it for all updates
c = colorbar;
caxis([Ttop Tbot]);

%*****  Solve Model Equations ******************************************

% Time-stepping loop
t = 0;  % initial time [years]
while t <= tend
    % Time increment (1 year per step)
    t = t + dt_yr;

    % Calculate the new temperature distribution (using diffusion model)
    T_new = T;
    for i = 2:size(T, 1)-1
        for j = 2:size(T, 2)-1
            % Simple 2D diffusion model (central difference)
            T_new(i, j) = T(i, j) + alpha * dt_yr * ...
                ((T(i+1, j) - 2*T(i, j) + T(i-1, j)) / h^2 + ...
                 (T(i, j+1) - 2*T(i, j) + T(i, j-1)) / h^2);
        end
    end
    T = T_new;

    % Update the image data with the new temperature values
    set(h_img, 'CData', T);
    title(['Temperature Distribution at t = ', num2str(t, '%.1f'), ' years']);
    drawnow;

    % Pause to allow dynamic update visualization
    pause(0.2);  % Pause for 0.2 seconds to show each year's update
end
