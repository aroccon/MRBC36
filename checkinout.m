clear all
nx=64;
ny=nx;
dx=2*pi/(nx-1);
dxi=1/dx;

fid = fopen('in.dat');
var = fread(fid, nx*ny, 'double');
fclose(fid);
in=reshape(var,[nx ny]);

fid = fopen('out.dat');
var = fread(fid, nx*ny, 'double');
fclose(fid);
out=reshape(var,[nx ny]);


figure(1)
subplot(1,2,1)
contourf(in)

subplot(1,2,2)
contourf(out)

