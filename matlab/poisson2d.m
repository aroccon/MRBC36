clear; close all; clc

% Parameters
nx = 128;
ny = nx;
Lx = 2*pi;
Ly = 1.0;
dx = Lx / nx;
dy = Ly / (ny - 1);

x = (0:nx-1) * dx;
y = linspace(0, Ly, ny);

n = 2;
m = 2;

% Allocate arrays in (x,y) layout
f        = zeros(nx, ny);   % f(i,j): x along rows, y along columns
u_exact  = zeros(nx, ny);

% Coefficient from analytical solution
A_coef = -1 / (n^2 + (2*pi*m/Ly)^2);

% Fill RHS and exact solution without meshgrid
for i = 1:nx
    for j = 1:ny
        f(i,j) = sin(n * x(i)) * cos(2*pi*m * y(j) / Ly);
        u_exact(i,j) = A_coef * sin(n * x(i)) * cos(2*pi*m * y(j) / Ly);
    end
end

% FFT in x-direction (along rows)
f_hat = fft(f, [], 1);  % size: (nx, ny)
kx = [0:nx/2 -nx/2+1:-1] * (2*pi / Lx);
kx2 = kx.^2;

% Build 2nd-order FD matrix with Neumann BCs in y
Ay = zeros(ny, ny);
for j = 2:ny-1
    Ay(j,j-1) = 1;
    Ay(j,j)   = -2;
    Ay(j,j+1) = 1;
end
Ay(1,1) = -2; Ay(1,2) = 2;           % Neumann BC at y=0
Ay(ny,ny) = -2; Ay(ny,ny-1) = 2;     % Neumann BC at y=Ly
Ay = Ay / dy^2;
I = eye(ny);

% Solve for each Fourier mode
u_hat = zeros(nx, ny);
for i = 1:nx
    L = Ay - kx2(i) * I;
    rhs = squeeze(f_hat(i,:))';  % size (ny, 1)

    if kx2(i) == 0
        L(1,:) = 0;
        L(1,1) = 1;
        rhs(1) = 0;
    end

    u_hat(i,:) = (L \ rhs)';
end

% Inverse FFT in x (along rows)
u = real(ifft(u_hat, [], 1));

% Compute L2 error (interior only)
err = u - u_exact;
L2 = sqrt(sum(err(:,2:ny-1).^2, 'all') / (nx * (ny-2)));
disp(['L2 error: ', num2str(L2)]);

% Plot results
figure
subplot(1,2,1)
contourf(x, y, u', 30, 'EdgeColor', 'none'); 
title('Numerical'); xlabel('x'); ylabel('y');

subplot(1,2,2)
contourf(x, y, u_exact', 30, 'EdgeColor', 'none'); 
title('Analytical'); xlabel('x'); ylabel('y');
