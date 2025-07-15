clear; clc; close all

% Parameters
nx = 64;
ny = nx;
Lx = 2*pi;
Ly = 1.0;
dx = Lx / nx;
dy = Ly / (ny - 1);
dyi=1/dy;
x = (0:nx-1)*dx;
y = linspace(0, Ly, ny);

n = 1;
m = 2;

% Exact solution and RHS
A_coef  = -1 / (n^2 + (2*pi*m/Ly)^2);
pext    = zeros(nx, ny);
rhsp    = zeros(nx, ny);
for i = 1:nx
    for j = 1:ny
        pext(i,j) =  A_coef*sin(n*x(i))*cos(2*pi*m*y(j)/Ly);
        rhsp(i,j) =  sin(n*x(i))*cos(2*pi*m*y(j)/Ly);
    end
end

% FFT in x (real-to-complex)
rhspc = fft(rhsp, [], 1);
rhspc = rhspc(1:(nx/2+1), :);  % Keep only positive frequencies

% Wavenumbers and their squares
kx = (0:nx/2) * (2*pi / Lx);
kx2 = kx.^2;

% Preallocate solution
pc = zeros(nx/2+1, ny);


for i = 1:nx/2+1
    % Set up tridiagonal system for each kx
    % The system is: a(j) * pc(i,j-1) + b(j) * pc(i,j) + c(j) * pc(i,j+1) = d(j)
    
    % Initialize diagonals and RHS
    a = zeros(ny,1);
    b = zeros(ny,1);
    c = zeros(ny,1);
    d = zeros(ny,1);
    
    for j = 1:ny
        a(j) =  1.0 * dyi^2;
        b(j) = -2.0 * dyi^2 - kx2(i);
        c(j) =  1.0 * dyi^2;
        d(j) = rhspc(i,j);
    end

    % Neumann BC at j = 1 (bottom)
    b(1) = -2.0 * dyi^2 - kx2(i);
    c(1) =  2.0 * dyi^2;
    a(1) =  0.0;


    % Neumann BC at j = ny (top)
    a(ny) =  2.0 * dyi^2;
    b(ny) = -2.0 * dyi^2 - kx2(i);
    c(ny) =  0.0;


       % Special handling for kx = 0 (mean mode)
    if (kx(i) == 0) 
        for j = 1:ny
            pc(i,j) = 0.0d0;
        end
    else
        % Thomas algorithm (TDMA)
        % Forward sweep
        for j = 2:ny
            factor = a(j) / b(j-1);
            b(j) = b(j) - factor * c(j-1);
            d(j) = d(j) - factor * d(j-1);
        end

        % Back substitution
        sol = zeros(ny,1);
        sol(ny) = d(ny) / b(ny);
        for j = ny-1:-1:1
            sol(j) = (d(j) - c(j) * sol(j+1)) / b(j);
        end

        % Store solution
        for j = 1:ny
            pc(i,j) = sol(j);
            %sol(j)
        end
    end
end

% Reconstruct full spectrum using Hermitian symmetry
p_hat_full = zeros(nx, ny);
p_hat_full(1:nx/2+1, :) = pc;

% Fill in the negative frequencies (complex conjugate symmetry)
for i = 2:nx/2
    p_hat_full(nx - i + 2, :) = conj(pc(i, :));
end


% Inverse FFT
p = real(ifft(p_hat_full, [], 1));


% L2 error ignoring boundaries
l2norm=0.d0;
for i=1:nx
    for j=1:ny	
        err = p(i,j)-pext(i,j);
        l2norm = l2norm + err*err;
    end
end
l2norm=sqrt(l2norm/nx/ny);
disp(['L2 error: ', num2str(l2norm)]);

figure(1)
clf
subplot(1,2,1)
hold on
surf(real(rhspc))
colorbar
hold off

subplot(1,2,2)
hold on
surf(imag(rhspc))
colorbar
hold off

% Plot results from Matlab
figure(3)
subplot(2,4,1)
contourf(x,y,rhsp', 30, 'EdgeColor', 'none');
hold on
title('RHS - Matlab'); 
xlabel('x'); 
ylabel('y');
colorbar
hold off

subplot(2,4,2)
contourf(x, y, p', 30, 'EdgeColor', 'none');
hold on
title('Numerical - Matlab'); 
xlabel('x'); 
ylabel('y');
colorbar
hold off

subplot(2,4,3)
contourf(x, y, pext', 30, 'EdgeColor', 'none');
title('Analytical'); 
xlabel('x'); 
ylabel('y');
colorbar
hold off

subplot(2,4,4)
contourf(x, y, (p-pext)', 30, 'EdgeColor', 'none');
title('Error - Matlab'); 
xlabel('x'); 
ylabel('y');
colorbar
hold off

