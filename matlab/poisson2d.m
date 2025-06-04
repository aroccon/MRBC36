% Parameters
Nx = 128;
Ny = 128;
Lx = 2 * pi;
Ly = 1.0;
dx = Lx/(Nx-1);
dy = Ly/(Ny - 1);
x = (0:Nx-1) * dx;
y = (0:Ny-1) * dy;
[X, Y] = meshgrid(x, y);

% Define RHS
f = sin(X) .* cos(pi * Y / Ly);

% Analytical solution
A = -1 / (1 + (pi / Ly)^2);
u_exact = A * sin(X) .* cos(pi * Y / Ly);

% FFT in x-direction
f_hat = fft(f, [], 2);

% Wavenumbers squared
kx = [0:Nx/2 -Nx/2+1:-1] * (2*pi/Lx);
kx2 = kx.^2;

% Build finite difference matrix with Neumann BCs (Ny x Ny)
A = zeros(Ny, Ny);
for j = 2:Ny-1
    A(j, j-1) = 1;
    A(j, j)   = -2;
    A(j, j+1) = 1;
end
A(1,1) = -2; 
A(1,2) = 2;
A(Ny,Ny-1) = 2; 
A(Ny,Ny) = -2;
A = A / dy^2;

% Identity matrix
I = eye(Ny);

% Solve for each Fourier mode
u_hat = zeros(Ny, Nx);
for i = 1:Nx
    % System matrix for this mode
    L = A - kx2(i) * I;

    % Right-hand side
    rhs = f_hat(:,i); 

    % Handle the singular kx = 0 mode separately
    if kx2(i) == 0
        % Set mean to zero or remove null space
        L(1,:) = 0;
        L(1,1) = 1;
        rhs(1) = 0;
    end

    % Gaussian elimination (basic implementation)
    % Forward elimination
    for k = 1:Ny-1
        if abs(L(k,k)) < 1e-14
            error('Pivot too small at row %d, matrix may be singular', k);
        end
        factor = L(k+1,k) / L(k,k);
        L(k+1,k:end) = L(k+1,k:end) - factor * L(k,k:end);
        rhs(k+1) = rhs(k+1) - factor * rhs(k);
    end

    % Back substitution
    u_col = zeros(Ny,1);
    u_col(Ny) = rhs(Ny) / L(Ny,Ny);
    for k = Ny-1:-1:1
        u_col(k) = (rhs(k) - L(k,k+1:end) * u_col(k+1:end)) / L(k,k);
    end

    % Store result
    u_hat(:,i) = u_col;
end

% Inverse FFT to get solution
u = real(ifft(u_hat, [], 2));

error = max(abs(u(:) - u_exact(:)), [], 'all');
disp(['Max error: ', num2str(error)]);


% Plot
figure(1)
subplot(1,2,1)
contourf(X, Y, u, 'EdgeColor', 'none');
xlabel('x'); 
ylabel('y'); 

subplot(1,2,2)
contourf(X, Y, u_exact, 'EdgeColor', 'none');
xlabel('x'); 
ylabel('y'); 