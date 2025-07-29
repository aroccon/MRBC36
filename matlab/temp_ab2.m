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
dt=0.06;
alpha=1.0;
beta=0.0;
tfin=1000;
ntmax=tfin/dt;

temp=zeros(nx,ny);
rhstemp=zeros(nx,ny);
rhstemp_o=zeros(nx,ny);

% theoretical for EE
dtmax=min(dx^2/4/diff,dy^2/4/diff);
disp(dtmax)
tic;

for t=1:ntmax
    for i=1:nx
        for j=2:ny-1
            ip=i+1;
            im=i-1;
            jp=j+1;
            jm=j-1;
            if (ip > nx) 
                ip=1;
            end
            if (im < 1) 
                im=nx;
            end
            rhstemp(i,j)= diff*((temp(ip,j)-2.d0*temp(i,j)+temp(im,j))*ddxi + (temp(i,jp) -2.d0*temp(i,j) +temp(i,jm))*ddyi);      
        end
    end

    for i=1:nx
        for j=2:ny-1
            % AB2
            temp(i,j) = temp(i,j)  + dt*(alpha*rhstemp(i,j)-beta*rhstemp_o(i,j));
            rhstemp_o(i,j)=rhstemp(i,j);
        end
    end
    alpha=1.0;
    beta=0.0;

    % impose BC on the temperature field
    for i=1:nx
        temp(i,1) =  1;
        temp(i,ny) = 0;
    end
end
elapsed = toc;
fprintf('Elapsed time: %.4f seconds\n', elapsed);

figure(1)
contourf(temp')