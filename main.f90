program main

use cufft
use iso_c_binding
use openacc 

implicit none
integer, parameter :: nx=64, ny=64
double precision, parameter :: pi=3.141592653589793d0
double precision :: dx, dy, lx, ly, acoeff
integer :: i, j, k, n, m
double precision, allocatable :: x(:), y(:)
double precision, allocatable :: rhsp(:,:),  p(:,:), pext(:,:)
double complex, allocatable :: rhspc(:,:)

! cufft plans
integer :: planf, planb, status


allocate(x(nx),y(ny))
allocate(rhsp(nx,ny),pext(nx,ny))
allocate(rhspc(nx/2+1,ny))

! Define domain size
lx = 2.d0*pi
ly = 1.d0
! Grid spacing: dx divided by nx to have perfect periodicity
dx = lx/nx 
dy = ly/(ny-1)
! Parameters for test solution
n= 1
m= 2
acoeff= -1.d0/((n*n + (2*pi*m/ly)**2))
! Compute axis
do i=1,nx-1
  x(i+1)=x(i) + dx
enddo
do j=1,ny
  y(j+1)=y(j) + dy
enddo

! Initialize rhsp term 
do j=1,ny
  do i=1,nx
    rhsp(i,j) = sin(n * x(i))*cos(2*pi*m * y(j)/ly)
    pext(i,j) = acoeff*rhsp(i,j)
  end do
end do

! create plans
! Forward along x
status=0
status=status+cufftCreate(planf)
status=status + cufftPlanMany(planf, 1, nx, [1], 1, nx, [1], 1, nx/2+1, CUFFT_D2Z, ny) 
if (status.ne.0) write(*,*) "Error in cuFFT plan FWD:", status

! Backward along x
status=0
status=status+cufftCreate(planb)
status=status + cufftPlanMany(planb, 1, nx, [1], 1, nx/2+1, [1], 1, nx, CUFFT_Z2D, ny) 

! write input field
open(unit=55,file='in.dat',form='unformatted',position='append',access='stream',status='new')
write(55) rhsp(:,:)
close(55)

! === Execute Forward FFT ===
status = cufftExecD2Z(planf, rhsp, rhspc)
if (status.ne.0) stop "cufftExecD2Z failed"

! === Execute Backward FFT ===
status = cufftExecZ2D(planb, rhspc, rhsp)
if (status.ne.0) stop "cufftExecZ2D failed"

rhsp=rhsp/nx


! write out field
open(unit=55,file='out.dat',form='unformatted',position='append',access='stream',status='new')
write(55) rhsp(:,:)
close(55)

deallocate(x,y)
deallocate(rhsp,p,pext,rhspc)

end program main
