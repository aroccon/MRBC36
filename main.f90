
program main

use cufft
use iso_c_binding
use openacc 

implicit none
integer, parameter :: nx=64, ny=64
double precision, parameter :: pi=3.141592653589793d0
double precision :: dx, dy, lx, ly, acoeff, q, l2norm, err, dyi, factor
integer :: i, j, k, n, m, t
double precision, allocatable :: x(:), y(:), kx(:), kx2(:)
double precision, allocatable :: u(:,:),  v(:,:), rhsu(:,:), rhsv(:,:)
double precision, allocatable :: rhsp(:,:),  p(:,:), pext(:,:)
double precision, allocatable :: rhsphi(:,:), phi(:,:), normx(:,:), normy(:,:)
double complex, allocatable :: rhspc(:,:), pc(:,:), rhs(:)
double complex :: a(ny), b(ny), c(ny), d(ny), sol(ny)
integer :: planf, planb, status, ntmax
double precision :: radius, eps, gamma, rho, mu, sigma

! declare parameters
ntmax=10
radius=0.5d0
eps=dx
gamma=0.0d0
rho=1.d0
mu=1.d0
sigma=0.d0

! allocate variables (avoid allocate/deallocate at run time)
allocate(x(nx),y(ny))
allocate(rhsp(nx,ny),p(nx,ny),pext(nx,ny))
allocate(rhspc(nx/2+1,ny),pc(nx/2+1,ny))
allocate(kx(nx/2+1))
allocate(kx2(nx/2+1))


! phase-field variables
allocate(phi(nx,ny),rhsphi(nx,ny))
allocate(normx(nx,ny),normy(nx,ny))




!##########################################################
! Define some basic quantities
!##########################################################
! Define domain size
lx = 2.d0*pi
ly = 1.d0
! Grid spacing: dx divided by nx to have perfect periodicity
dx = lx/nx 
dy = ly/(ny-1)
dyi = 1.d0/dy
dxi=1.d0/dx
ddxi=dxi*dxi
ddyi=dyi*dyi

! Parameters for test solution (debug only)
!n= 1
!m= 2
!acoeff= -1.d0/((n*n + (2*pi*m/ly)**2))
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

! write input field (for debug)
!open(unit=55,file='in.dat',form='unformatted',position='append',access='stream',status='new')
!write(55) rhsp(:,:)
!close(55)


!##########################################################
! Init fields
!##########################################################
do i=1,nx
  do j=1,ny
    phi(i,j)=1.0d0
  enddo
enddo


!##########################################################
! Start temporal loop
!##########################################################


!do t=1,ntmax
!##########################################################
! Advance phase-field (constant ATM)
!##########################################################

  ! advection
  do i=1,nx
    do j=2,ny-1
      ip=i+1
      im=im-1
      jp=j+1
      jm=j-1
      if (ip .gt. nx) ip=1
      if (im .lt. 1) im=nx
      rhsphi(i,j)=-(u(ip,j)*0.5*(phi(ip,j)+phi(i,j)) - u(i,j)*0.5*(phi(i,j)+phi(im,j)))*dxi -(v(i,jp)*0.5*(phi(i,jp)+phi(i,j)) - v(i,j)*0.5*(phi(i,j)+phi(i,jm)))*dyi;
    enddo
  enddo
  ! diffusion
  do i=1:nx
    do j=2:ny-1
      ip=i+1
      im=im-1
      jp=j+1
      jm=j-1
      if (ip .gt. nx) ip=1
      if (im .lt. 1) im=nx
      rhsphi(i,j)=rhsphi(i,j) + gams*(eps*(phi(ip,j)-2*phi(i,j)+phi(im,j))*ddxi +eps*(phi(i,jp) -2*phi(i,j) +phi(i,jm))*ddyi);      
    enddo
  enddo
  !compute normals
  do i=1,nx
    do j=2,ny-1
      ip=i+1
      im=im-1
      jp=j+1
      jm=j-1
      if (ip .gt. nx) ip=1
      if (im .lt. 1) im=nx
      normx(i,j) = (phi(ip,j)-phi(im,j));
      normy(i,j) = (phi(i,jp)-phi(i,jm));
      normod = 1.0/(sqrt(normx(i,j)^2+normy(i,j)^2) + 1.e-16);
      normx(i,j) = normx(i,j)*normod;
      normy(i,j) = normy(i,j)*normod;
    enddo
  enddo
  !add sharpening
  do i=1,nx
    do j=2,ny-1
      ip=i+1
      im=im-1
      jp=j+1
      jm=j-1
      if (ip .gt. nx) ip=1
      if (im .lt. 1) im=nx
      rhsphi(i,j)=rhsphi(i,j)+gams*(((phi(ip,j)^2-phi(ip,j))*normx(ip,j) - (phi(im,j)^2-phi(im,j))*normx(im,j))*0.5*dxi + ((phi(i,jp)^2-phi(i,jp))*normy(i,jp) - (phi(i,jm)^2-phi(i,jm))*normy(i,jm))*0.5*dyi);
    enddo
  enddo
  !%boiling term
  !for i=1:nx
  !   for j=1:ny
  !       melt(i,j) = phi(i,j)*(1-phi(i,j))/eps*mb;
  !       rhsphi(i,j)= rhsphi(i,j) + (2-rhor)*melt(i,j)*((1-phi(i,j))/rhol + phi(i,j)/rhov);
  !   end
  !end
  ! phase-field n+1
  do i=1,nx
    do j=2,ny-2
      phi(i,j) = phi(i,j) + dt*rhsphi(i,j);
    enddo
  enddo
  ! impose BC on the phase-field
  do i=1,nx
    phi(i,1) = phi(i,2)
    phi(i,ny) = phi(i,ny-1)
  enddo
  ! no flux at the walls


!##########################################################
!Start of Poisson solver, pressure in physical space obtained
!##########################################################

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
    ! Thomas algorithm (TDMA) for tridiagonal system 
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

! write pressure (debug only)
open(unit=55,file='out.dat',form='unformatted',position='append',access='stream',status='new')
write(55) p(:,:)
close(55)
!##########################################################
!End of Poisson solver, pressure in physical space obtained
!##########################################################


!##########################################################
!Start correction step
!##########################################################


deallocate(x,y)
deallocate(a,b,c,d,sol)
deallocate(kx,kx2)
deallocate(rhsp,p,rhspc)

end program main
