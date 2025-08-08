%clear all
lx=2.0;
ly=1.0;
nx = 512; 
ny = 256;
dx = lx/nx; 
dy = ly/(ny-1);
dt = 0.01;
kappa = 0.1;
ntmax= 99;

r_y = kappa*dt/dy^2;
r_x = kappa*dt/dx^2;

% Initialize temperature T and velocity fields u,v (example)
T = zeros(nx, ny);
T(:,1) = 0.0;  % top cold
T(:,ny) = 1.0; % bottom hot

u = zeros(nx, ny); % example constant velocity in x
v = zeros(nx, ny);    % zero velocity in y

% Preallocate
Adv = zeros(nx, ny);
RHS = zeros(nx, ny);
T_new = zeros(nx, ny);

for t=1:ntmax
    % Compute explicit advection term Adv(i,j)
    for i = 1:nx
        ip = mod(i, nx) + 1;   % periodic x+
        im = mod(i-2, nx) + 1; % periodic x-
        for j = 2:ny-1
            jp = j + 1;
            jm = j - 1;
            adv_x = -(u(ip,j)*0.5*(T(ip,j)+T(i,j)) - u(i,j)*0.5*(T(i,j)+T(im,j))) / dx;
            adv_y = -(v(i,jp)*0.5*(T(i,jp)+T(i,j)) - v(i,j)*0.5*(T(i,j)+T(i,jm))) / dy;
            Adv(i,j) = adv_x + adv_y;
        end
    end

    % Build RHS for CN implicit y diffusion solve
    for i = 1:nx
        ip = mod(i, nx) + 1;   % periodic x+
        im = mod(i-2, nx) + 1; % periodic x-
    
        for j = 2:ny-1
            % Diffusion in x explicit
            diff_x_explicit = kappa * dt * (T(ip,j) - 2*T(i,j) + T(im,j)) / dx^2;
            % Add half of old y-diffusion (Crank-Nicolson)
            diff_y_old = kappa * (T(i,j+1) - 2*T(i,j) + T(i,j-1))/dy/dy;
            RHS(i,j) = T(i,j) + dt * (-Adv(i,j) + diff_x_explicit) + 0.5*dt*diff_y_old;
        end
        % Enforce BC on RHS
        RHS(i,1) = 0.0;   % top cold
        RHS(i,ny) = 1.0;  % bottom hot
    end

    % Implicit solve in y-direction (TDMA) for each column i
    for i = 1:nx
        % TDMA coefficients for Crank-Nicolson
        a = -0.5*r_y * ones(ny,1); % sub-diagonal
        b = (1 + r_y) * ones(ny,1); % diagonal
        c = -0.5*r_y * ones(ny,1); % super-diagonal
    
        % Adjust for Dirichlet BC rows
        b(1) = 1.0; c(1) = 0.0;
        a(ny) = 0.0; b(ny) = 1.0;
    
        % TDMA forward sweep
        cp = zeros(ny,1);
        dp = zeros(ny,1);
        cp(1) = c(1)/b(1);
        dp(1) = RHS(i,1)/b(1);
    
        for j = 2:ny
            denom = b(j) - a(j)*cp(j-1);
            cp(j) = c(j)/denom;
            dp(j) = (RHS(i,j) - a(j)*dp(j-1))/denom;
        end
    
        % TDMA backward substitution
        T_new(i,ny) = dp(ny);
        for j = ny-1:-1:1
            T_new(i,j) = dp(j) - cp(j)*T_new(i,j+1);
        end
    end

    % Update temperature field
    T = T_new;
end

% Plot updated temperature
imagesc(T');
colorbar
title('Final temperature')
axis equal tight
