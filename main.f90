
program main

use cufft
use iso_c_binding
use openacc 

implicit none
integer, parameter :: nx=128, ny=128
double precision, parameter :: pi=3.141592653589793d0
double precision :: dx, dy, lx, ly, acoeff, q, l2norm, err, dyi, factor
integer :: i, j, k, n, m
double precision, allocatable :: x(:), y(:), kx(:), kx2(:)
double precision, allocatable :: rhsp(:,:),  p(:,:), pext(:,:)
double complex, allocatable :: rhspc(:,:), pc(:,:), rhs(:)
double complex :: a(ny), b(ny), c(ny), d(ny), sol(ny)

! cufft plans
integer :: planf, planb, status


allocate(x(nx),y(ny))
allocate(rhsp(nx,ny),p(nx,ny),pext(nx,ny))
allocate(rhspc(nx/2+1,ny),pc(nx/2+1,ny))
allocate(kx(nx/2+1))
allocate(kx2(nx/2+1))

! Define domain size
lx = 2.d0*pi
ly = 1.d0
! Grid spacing: dx divided by nx to have perfect periodicity
dx = lx/nx 
dy = ly/(ny-1)
dyi = 1.d0/dy
! Parameters for test solution
n= 1
m= 2
acoeff= -1.d0/((n*n + (2*pi*m/ly)**2))
! Compute axis
x(1)=0.d0
y(1)=0.d0
do i=1,nx-1
  x(i+1)=x(i) + dx
enddo
do j=1,ny-1
  y(j+1)=y(j) + dy
enddo

! --- Wavenumbers and squares ---
do i = 1,nx/2+1
  kx(i) = (i-1)*(2*pi/Lx)
  kx2(i) = kx(i)**2
end do

! Initialize rhsp term 
do j=1,ny
  do i=1,nx
    rhsp(i,j) = sin(n * x(i))*cos(2*pi*m*y(j)/ly)
    pext(i,j) = acoeff*rhsp(i,j)
  end do
end do

!write(*,*) "x", x
!write(*,*) "y", y

! create plans
! Forward along x
status=0
status=status + cufftCreate(planf)
status=status + cufftPlanMany(planf, 1, nx, [nx, ny], 1, nx, [nx/2+1,ny], 1, nx/2+1, CUFFT_D2Z, ny) 
if (status.ne.0) write(*,*) "Error in cuFFT plan FWD:", status

! Backward along x
status=0
status=status + cufftCreate(planb)
status=status + cufftPlanMany(planb, 1, nx, [nx/2+1,ny], 1, nx/2+1, [nx,ny], 1, nx, CUFFT_Z2D, ny) 

! write input field
open(unit=55,file='in.dat',form='unformatted',position='append',access='stream',status='new')
write(55) rhsp(:,:)
close(55)

! === Execute Forward FFT ===
!$acc data copy(rhsp,rhspc)
!$acc host_data use_device(rhsp,rhspc)
status = cufftExecD2Z(planf, rhsp, rhspc)
!$acc end host_data
!$acc end data
if (status.ne.0) stop "cufftExecD2Z failed"

!$acc kernels
do i=1,nx/2+1
  ! Set up tridiagonal system for each kx
  ! The system is: (A_j) * pc(kx,j-1) + (B_j) * pc(kx,j) + (C_j) * pc(kx,j+1) = rhs(kx,j)
  ! FD2: -pc(j-1) + 2*pc(j) - pc(j+1)  --> Laplacian in y
  ! Neumann BC: d/dy pc = 0 at j=1 and j=ny

  ! Fill diagonals and rhs for each y
  do j = 1, ny
    a(j) =  1.0d0*dyi*dyi
    b(j) = -2.0d0*dyi*dyi - kx2(i)
    c(j) =  1.0d0*dyi*dyi
    d(j) =  rhspc(i,j)
  end do

  ! Neumann BC at j=1 (bottom)
  b(1) = -2.0d0*dyi*dyi - kx2(i)
  c(1) =  2.0d0*dyi*dyi
  a(1) =  0.0d0

  ! Neumann BC at j=ny (top)
  a(ny) =  2.0d0*dyi*dyi
  b(ny) = -2.0d0*dyi*dyi - kx2(i)
  c(ny) =  0.0d0

  ! Special handling for kx=0 (mean mode)
  if (kx(i) == 0.0d0) then
    do j = 1, ny
        pc(i,j) = 0.0d0
    end do
  else
    ! Thomas algorithm (TDMA) for tridiagonal system [1][5]
    ! Forward sweep
    do j = 2, ny
      factor = a(j)/b(j-1)
      b(j) = b(j) - factor * c(j-1)
      d(j) = d(j) - factor * d(j-1)
    end do

    ! Back substitution
    sol(ny) = d(ny) / b(ny)
    do j = ny-1, 1, -1
      sol(j) = (d(j) - c(j) * sol(j+1)) / b(j)
    end do

    ! Store solution
    do j = 1, ny
      pc(i,j) = sol(j)
      !write(*,*) "sol", sol(j)
    end do
  end if
end do
!$acc end kernels



! === Execute Backward FFT ===
!$acc data copy(pc,p)
!$acc host_data use_device(pc,p)
status = cufftExecZ2D(planb, pc, p)
!$acc end host_data
!$acc end data
if (status.ne.0) stop "cufftExecZ2D failed"

!$acc kernels
p=p/nx
!$acc end kernels

!write(*,*) "max p", maxval(real(pc))

! write out analytical solution
!open(unit=55,file='pext.dat',form='unformatted',position='append',access='stream',status='new')
!write(55) pext(:,:)
!close(55)

! write out field
open(unit=55,file='out.dat',form='unformatted',position='append',access='stream',status='new')
write(55) p(:,:)
close(55)



l2norm=0.d0
do i=1,nx
 do j=1,ny	
  err = p(i,j)-pext(i,j)
  l2norm = l2norm + err*err
 enddo
enddo
	
l2norm=sqrt(l2norm/nx/ny)

write(*,*) "l2norm", l2norm

deallocate(x,y)
!deallocate(a,b,c,d,sol)
!deallocate(kx,kx2)
!deallocate(rhsp,p,pext,rhspc)

end program main
