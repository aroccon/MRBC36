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


figure(2)
clf
subplot(1,2,1)
hold on
surf(out)
colorbar
hold off

subplot(1,2,2)
hold on
surf(in)
colorbar
hold off


