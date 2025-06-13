program main

use cufft
use iso_c_binding
use openacc 

implicit none
integer, parameter :: nx=64, ny=64
double precision, parameter :: pi=3.141592653589793d0
double precision :: dx, dy, lx, ly, acoeff
integer :: i, j, k, n, m
double precision, allocatable :: x(:), y(:), kx(:), kx2(:)
double precision, allocatable :: rhsp(:,:),  p(:,:), pext(:,:)
double complex, allocatable :: rhspc(:,:)
double precision, alloctable :: a(:), b(:), c(:), d(:), sol(:)

! cufft plans
integer :: planf, planb, status


allocate(x(nx),y(ny))
allocate(rhsp(nx,ny),p(nx,ny),pext(nx,ny))
allocate(rhspc(nx/2+1,ny))
allocate(kx(nx/2+1))
allocate(kx2(nx/2+1))
allocate(a(ny),b(ny),c(ny),d(ny),sol(ny))

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


! Arrays a, b, c, d represent the sub-diagonal, diagonal, super-diagonal, and RHS vectors respectively.
! Solve for each Fourier mode kx
do i=1,nx/2+1
  kx = 2.0*pi*(k-1)/Lx

  ! Setup tridiagonal system for y
  do j=1,ny
    a(j) = -1.0
    b(j) = 2.0 + kx2(i)*dy*dy
    c(j) = -1.0
    d(j) = dy*dy*real(rhspc(i,j))
  end do

  ! Modify for Neumann BC at boundaries
  b(1) = 1.0
  c(1) = -1.0
  d(1) = 0.0

  a(ny) = -1.0
  b(ny) = 1.0
  d(ny) = 0.0

  ! ==TDMA forward elimination ---
  do j=2,ny
    q = a(j)/b(j-1)
    b(j) = b(j) - q*c(j-1)
    d(j) = d(j) - q*d(j-1)
  end do

  ! == TDMA Back substitution ---
  sol(ny) = d(ny)/b(ny)
  do j=ny-1,1,-1
    sol(j) = (d(j) - c(j)*sol(j+1))/b(j)
  end do

  ! Store solution in spectral space
  do j=1,ny
    pc(i,j) = cmplx(sol(j), 0.0)
  end do
end do

! === Execute Backward FFT ===
status = cufftExecZ2D(planb, pc, p)
if (status.ne.0) stop "cufftExecZ2D failed"


! write out field
open(unit=55,file='out.dat',form='unformatted',position='append',access='stream',status='new')
write(55) rhsp(:,:)
close(55)

deallocate(x,y)
deallocate(a,b,c,d,sol)
deallocate(kx,kx2)
deallocate(rhsp,p,pext,rhspc)

end program main
