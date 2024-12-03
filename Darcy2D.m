%*****  1D ADVECTION DIFFUSION MODEL OF HEAT TRANSPORT  *******************

%*****  Initialise Model Setup

% create x-coordinate vectors
xc = h/2:h:W-h/2;      % x-coordinate vector for cell centre positions [m]
zc = h/2:h:D-h/2;      % z-coordinate vector for cell centre positions [m]
xf = 0:h:W;            % x-coordinate vectore for cell face positions [m]
zf = 0:h:D;            % z-coordinate vectore for cell face positions [m]
[Xc,Zc] = meshgrid(xc,zc);  % create 2D coordinate arrays

% set up index array for boundary conditions
ix3 = [          Nx,1:Nx,1       ];  % periodic sides
ix5 = [     Nx-1,Nx,1:Nx,1 ,2    ];
ix7 = [Nx-2,Nx-1,Nx,1:Nx,1 ,2 ,3 ];
iz3 = [           1,1:Nz,Nz      ];  % closed/insulating top/bot
iz5 = [        1, 1,1:Nz,Nz,Nz   ];
iz7 = [   1,   1, 1,1:Nz,Nz,Nz,Nz];

% create smooth random perturbation field
rng(15);
dr  = randn(Nz,Nx);
for ii = 1:10
    dr  = dr + (diff(dr(iz3,:),2,1) + diff(dr(:,ix3),2,2))/8;
end

% set initial condition for temperature at cell centres
T   = Ttop + (Tbot-Ttop)./D.*Zc + dr*1;  % initialise T array on linear gradient

% initialise density and mobility
rho = rho0.*(1 - aT.*(T-Ttop));
KD  = KD0.*ones(Nz,Nx);

% initialise Darcy velocity and non-hydrostatic pressure
u = zeros(Nz,Nx+1);
w = zeros(Nz+1,Nx);
p = zeros(Nz,Nx);

% initialise residual and update fields
res_T = 0*T;  upd_T = 0*T;  dTdt = 0*T;
res_p = 0*p;  upd_p = 0*p;

% initialise output figure with initial condition
figure(1); clf
makefig(xc,zc,T,p,w,u,0,yr)

%*****  Solve Model Equations

dt = CFL * min((h/2)/max([w(:);u(:)]),(h/2)^2/max(kT(:))); % initial time step [s]
t  = 0;  % initial time [s]
k  = 0;  % initial time step count

% loop through time steps until stopping time reached
while t <= tend

    % increment time and step count
    t = t+dt;
    k = k+1;

    % print time step header
    fprintf(1,'\n\n*****  step = %d;  dt = %1.3e;  time = %1.3e \n\n',k,dt,t)

    % store old temperature and rate
    dTdto = dTdt;
    To    = T;

    % get time step step size
    dt    = CFL * min((h/2)/max([w(:);u(:)]),(h/2)^2/max(kT(:))); % time step [s]

    resnorm = 1;  % initialise residual norm
    it      = 0;  % initialise iteration count

    % loop through pseudo-transient iterations until convergence criterion reached
    while resnorm > tol

        % update temperature every 'nup' iterations
        if ~mod(it,nup) && k>1

            % get rate of change
            dTdt = diffusion(T,kT,h,ix3,iz3) + advection(T,u,w,h,ix7,iz7,ADVN);

            % get temperature residual
            res_T = (T - To)/dt - (dTdt + dTdto)/2;

            % set isothermal boundaries on top/bot
            res_T(1  ,:) = 0;
            res_T(end,:) = 0;

            % get solution update
            upd_T = - alpha*res_T*dt/2;

            % update solution
            T     = T + upd_T;

        end

        % get density and density contrast
        rho   = rho0.*(1 - aT.*(T-Ttop));  % T-dependent density
        Drho  = rho - mean(rho,2);         % subtract horizontal mean
        Drhoz = (Drho(iz3(1:end-1),:)+Drho(iz3(2:end),:));  % on z-faces
        Drhoz([1 end],:) = 0;              % no flow across top/bot bounds

        % get p-dependent segregation mobility
        KD  = KD0 + c.*max(0,p).^m;   % nonlinear form
        KDx = 0.5 * (KD(:, ix5(1:Nx)) + KD(:, ix5(2:Nx+1)));  % on x-faces
        KDz = 0.5 * (KD(iz7(1:Nz), :) + KD(iz7(2:Nz+1), :));  % on z-faces

        % get Darcy flux (vD = - KD (Grad(p) - D(rho) g)
        u = - KDx .* (p(:, ix5(2:Nx+1)) - p(:, ix5(1:Nx))) / h - D(rho) * g0;  % x-speed
        w = - KDz .* (p(iz7(2:Nz+1), :) - p(iz7(1:Nz), :)) / h - D(rho) * g0;  % z-speed

        % get pressure residual (Div.v = 0)
        res_p = ( (u(:, ix5(2:Nx+1)) - u(:, ix5(1:Nx))) / h + (w(iz7(2:Nz+1), :) - w(iz7(1:Nz), :)) / h );

        % get pseudo-time step size
        max_u_w = max(abs(u(:)), abs(w(:)));  % 计算u和w的最大值
        dtau = min((h/2) / max_u_w, (h/2)^2 / max(kD(:)));

        % get solution update
        upd_p = - res_p * dtau;

        % update solution
        p = p + upd_p;

        it = it+1; % increment iteration count

        % get residual norm and print convergence every 'nup' iterations
        if ~mod(it,nup)
            resnorm = norm(upd_T(:),2)./norm(T(:)+eps,2) ...
                    + norm(upd_p(:),2)./norm(p(:)+eps,2);
            if isnan(resnorm); error('!!! Solver failed with nan !!!'); end
            fprintf(1,'     it = %d;  res = %e \n',it,resnorm); 
        end

    end

    % plot model progress every 'nop' time steps
    if ~mod(k,nop)
        makefig(xc,zc,T,p,w,u,t,yr);
    end

end


%*****  Utility Functions  ************************************************

% Function to make output figure
function makefig(x,z,T,p,w,u,t,yr)

clf; 

% plot temperature in subplot 1
subplot(2,2,1);
imagesc(x,z,T); axis equal tight; colorbar; hold on
contour(x,z,T,[100,150,200],'k');

ylabel('z [m]','FontSize',15)
title('Temperature [C]','FontSize',17)

% plot pressure in subplot 2
subplot(2,2,2)
imagesc(x,z,p); axis equal tight; colorbar

xlabel('x [m]','FontSize',15)
ylabel('z [m]','FontSize',15)
title('Pressure [Pa]','FontSize',17)

% plot z-speed in subplot 3
subplot(2,2,3)
imagesc(x,z,-w*yr); axis equal tight; colorbar

xlabel('x [m]','FontSize',15)
ylabel('z [m]','FontSize',15)
title('z-Speed [m/yr]','FontSize',17)

% plot x-speed in subplot 1
subplot(2,2,4)
imagesc(x,z,u*yr); axis equal tight; colorbar

xlabel('x [m]','FontSize',15)
ylabel('z [m]','FontSize',15)
title('x-Speed [m/yr]','FontSize',17)

sgtitle(['time = ',num2str(t/yr),' [yr]'],'FontSize',17)
drawnow;

end

% Function to calculate diffusion rate
function [dTdt] = diffusion(f, k, h, ix, iz)
    % Calculate the diffusion rate using central difference method
    % f: the field (e.g., temperature)
    % k: the diffusion coefficient
    % h: grid spacing
    % ix: index for x-boundary conditions
    % iz: index for z-boundary conditions

    % Initialize diffusion rate
    dTdt = zeros(size(f));

    % Central difference for diffusion (2nd order spatial derivative)
    % in x-direction (horizontal diffusion)
    dTdt(:, 2:end-1) = k .* (f(:, 3:end) - 2*f(:, 2:end-1) + f(:, 1:end-2)) / (h^2);

    % in z-direction (vertical diffusion)
    dTdt(2:end-1, :) = dTdt(2:end-1, :) + k .* (f(3:end, :) - 2*f(2:end-1, :) + f(1:end-2, :)) / (h^2);

    % Apply boundary conditions based on ix, iz (if applicable)
    dTdt(iz) = 0;  % e.g., zero-gradient boundary at top and bottom (insulating)
    dTdt(:, ix) = 0;  % e.g., periodic or insulating boundary at sides (adjust as needed)
    
    % Adjust for any other specific boundary conditions
end



% Function to calculate advection rate
function dTdt = advection(f, u, w, h, ix, iz, ADVN)
    % Calculate the advection rate using a first-order upwind method or another method like central difference
    % f: the field (e.g., temperature)
    % u: velocity in the x-direction
    % w: velocity in the z-direction
    % h: grid spacing
    % ix: index for x-boundary conditions
    % iz: index for z-boundary conditions
    % ADVN: a flag for which advection scheme to use (e.g., upwind or central)

    % Initialize advection rate
    dTdt = zeros(size(f));
    
    if ADVN == 1
        % Upwind scheme for advection (first-order accurate)
        
        % In x-direction (upwind method)
        dTdt(:, 2:end-1) = -u(:, 2:end-1) .* (f(:, 2:end-1) - f(:, 1:end-2)) / h;

        % In z-direction (upwind method)
        dTdt(2:end-1, :) = dTdt(2:end-1, :) - w(2:end-1, :) .* (f(2:end-1, :) - f(1:end-2, :)) / h;
        
    elseif ADVN == 2
        % Central difference method for advection (second-order accurate)

        % In x-direction
        dTdt(:, 2:end-1) = -u(:, 2:end-1) .* (f(:, 3:end) - f(:, 1:end-2)) / (2*h);

        % In z-direction
        dTdt(2:end-1, :) = dTdt(2:end-1, :) - w(2:end-1, :) .* (f(3:end, :) - f(1:end-2, :)) / (2*h);
        
    else
        error('Invalid advection scheme specified.');
    end

    % Apply boundary conditions based on ix, iz
    dTdt(iz) = 0;  % e.g., zero-gradient boundary at top and bottom (insulating)
    dTdt(:, ix) = 0;  % e.g., periodic or insulating boundary at sides (adjust as needed)

    % Adjust for other specific boundary conditions (e.g., fixed value or no flux)
end



% Function to calculate 5th order WENO polynomials from Jiang & Shu, 1996, J Comp Physics
function [fface] = weno5poly (fmm, fm, fc, fp, fpp)

% 5th order polynomials
p1 = (2*fmm - 7*fm + 11*fc )/6;
p2 = ( -fm  + 5*fc +  2*fp )/6;
p3 = (2*fc  + 5*fp -    fpp)/6;

% smoothness measure
b1 = 13/12*(fmm - 2*fm + fc ).^2 + 1/4*(  fmm - 4*fm + 3*fc ).^2;
b2 = 13/12*(fm  - 2*fc + fp ).^2 + 1/4*(  fm  -          fp ).^2;
b3 = 13/12*(fc  - 2*fp + fpp).^2 + 1/4*(3*fc  - 4*fp +   fpp).^2;

% weights
g   = [1/10, 6/10, 3/10];
wp1 = g(1)./(b1.^2 + eps);
wp2 = g(2)./(b2.^2 + eps);
wp3 = g(3)./(b3.^2 + eps);

% get flux (normalize weights at the same time)
fface = (wp1.*p1 + wp2.*p2 + wp3.*p3) ./ (wp1 + wp2 + wp3) ;

end