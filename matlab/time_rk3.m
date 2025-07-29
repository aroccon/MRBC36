% clear all
% Temperature advancement using RK3
nx=64;
ny=64;
lx=1;
ly=1;
dx=lx/nx;
dy=ly/ny;
ddxi=1/dx/dx;
ddyi=1/dy/dy;
diff=0.001;
dt=0.22;
tfin=1000;
ntmax=tfin/dt;

temp=zeros(nx,ny);
rhstemp=zeros(nx,ny);
rhs_old=zeros(nx,ny);

% theoretical for EE
dtmax=min(dx^2/4/diff,dy^2/4/diff);
disp(dtmax)
tic;

% RK3 coefficients (Wray low-storage)
alpha = [8/15, 5/12, 3/4];
beta  = [0, -17/60, -5/12];

% theoretical for EE
dtmax=min(dx^2/4/diff,dy^2/4/diff);
disp(dtmax)

tic;

for t=1:ntmax

    rhs_old(:,:) = 0;
    for stage = 1:3
        % ---- compute RHS ----
        for i=1:nx
            for j=2:ny-1
                ip=i+1; im=i-1; jp=j+1; jm=j-1;
                if (ip > nx), ip=1; end
                if (im < 1),  im=nx; end
                rhstemp(i,j)= diff*((temp(ip,j)-2*temp(i,j)+temp(im,j))*ddxi + ...
                                (temp(i,jp) -2*temp(i,j) +temp(i,jm))*ddyi);
            end
        end

        % update solution with alpha and beta weights on RHS
        for i=1:nx
            for j=2:ny-1
                temp(i,j) = temp(i,j) + dt * (alpha(stage)*rhstemp(i,j) + beta(stage)*rhs_old(i,j));
            end
        end

        % ---- enforce BC ----
        for i=1:nx
            temp(i,1) = 1;
            temp(i,ny) = 0;
        end
    end

    % store current RHS for next stage
    rhs_old = rhstemp;
end


elapsed = toc;
fprintf('Elapsed time: %.4f seconds\n', elapsed);

figure(1)
contourf(temp')