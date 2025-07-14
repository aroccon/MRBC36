
program main

use cufft
use iso_c_binding
use openacc 

implicit none
integer, parameter :: nx=64, ny=64
double precision, parameter :: pi=3.141592653589793d0
double precision :: dx, dy, lx, ly, acoeff, q, l2norm, err, dyi, factor
integer :: i, j, k, n, m, t
integer :: ip,im,jp,jm
double precision, allocatable :: x(:), y(:), kx(:), kx2(:)
double precision, allocatable :: u(:,:),  v(:,:), rhsu(:,:), rhsv(:,:)
double precision, allocatable :: rhsp(:,:),  p(:,:), pext(:,:)
double precision, allocatable :: rhsphi(:,:), phi(:,:), normx(:,:), normy(:,:)
double precision, allocatable :: rhstemp(:,:), temp(:,:)
double complex, allocatable :: rhspc(:,:), pc(:,:), rhs(:)
double complex :: a(ny), b(ny), c(ny), d(ny), sol(ny)
integer :: planf, planb, status, ntmax
double precision :: radius, eps, epsi, gamma, rho, mu, sigma, dxi, ddxi, ddyi, normod, dt, umax
double precision :: pos, epsratio, times, timef, difftemp

!##########################################################
! declare parameters
! Define some basic quantities
ntmax=10
radius=0.5d0
epsratio=1.0d0
gamma=0.0d0
rho=1.d0
mu=1.d0
sigma=0.d0
radius=0.3d0
difftemp=0.01d0
! Define domain size
lx = 1.d0
ly = 1.d0
! Grid spacing: dx divided by nx to have perfect periodicity
dx = lx/nx 
dy = ly/(ny-1)
dxi=1.d0/dx
dyi=1.d0/dy
ddxi=dxi*dxi
ddyi=dyi*dyi
dt=0.0001
eps=max(dx,dy)
epsi=1.d0/eps

!########################################################
! allocate variables (avoid allocate/deallocate at run time)
allocate(x(nx),y(ny))
allocate(rhsp(nx,ny),p(nx,ny),pext(nx,ny))
allocate(rhspc(nx/2+1,ny),pc(nx/2+1,ny))
allocate(kx(nx/2+1))
allocate(kx2(nx/2+1))

! velocity variables (defined on cell faces)
allocate(u(nx,ny),v(nx,ny-1))

! phase-field variables (defined on centers)
allocate(phi(nx,ny),rhsphi(nx,ny))
allocate(normx(nx,ny),normy(nx,ny))

! temperature variables (defined on centers)
allocate(temp(nx,ny),rhstemp(nx,ny))


!########################################################

write(*,*) "NEMESI36"
write(*,*) "Code for 2D phase-field simulation"
write(*,*) "Grid.     :", nx, ny
write(*,*) "Lx and Ly :", lx, ly
write(*,*) "Dx and Dy :", dx, dy
write(*,*) "Time step :", dt
write(*,*) "eps       :", eps

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
write(*,*) "Initialize velocity and phase-field"
! u velocity
do i=1,nx
  do j=1,ny
    u(i,j)=0.0d0
  enddo
enddo
! v velocity
do i=1,nx
  do j=1,ny-1
    v(i,j)=0.0d0
  enddo
enddo
! phase-field
do i=1,nx
  do j=1,ny
    pos=(x(i)-lx/2)**2d0 + (y(j)-ly/2)**2d0
    phi(i,j)=0.5d0*(1.d0-tanh((sqrt(pos)-radius)/2.d0/eps))
  enddo
enddo


!##########################################################
! Start temporal loop
!##########################################################

write(*,*) "Start temporal loop"
do t=1,ntmax
  call cpu_time(times)
  write(*,*) "Time step",t,"of",ntmax
  !##########################################################
  ! STEP 1: Advance phase-field 
  !##########################################################
  ! Advection + diffusion term
  umax=0.0d0 ! replace then with real vel max
  gamma=1.0d0*umax
  !$acc kernels
  do i=1,nx
    do j=2,ny-1
      ip=i+1
      im=im-1
      jp=j+1
      jm=j-1
      if (ip .gt. nx) ip=1
      if (im .lt. 1) im=nx
      ! Advection 
      rhsphi(i,j)=-(u(ip,j)*0.5*(phi(ip,j)+phi(i,j)) - u(i,j)*0.5*(phi(i,j)+phi(im,j)))*dxi -(v(i,jp)*0.5*(phi(i,jp)+phi(i,j)) - v(i,j)*0.5*(phi(i,j)+phi(i,jm)))*dyi
      ! Add diffusion
      rhsphi(i,j)=rhsphi(i,j) + gamma*eps*((phi(ip,j)-2.d0*phi(i,j)+phi(im,j))*ddxi + (phi(i,jp) -2.d0*phi(i,j) +phi(i,jm))*ddyi)     
      ! Compute normals
      normx(i,j) = (phi(ip,j)-phi(im,j))
      normy(i,j) = (phi(i,jp)-phi(i,jm))
      normod = 1.0d0/(sqrt(normx(i,j)**2d0 + normy(i,j)**2d0) + 1.e-16)
      normx(i,j) = normx(i,j)*normod
      normy(i,j) = normy(i,j)*normod 
    enddo
  enddo
  !$acc end kernels

  !$acc kernels
  ! compute sharpening (then move to ACDI)
  do i=1,nx
    do j=2,ny-1
      ip=i+1
      im=im-1
      jp=j+1
      jm=j-1
      if (ip .gt. nx) ip=1
      if (im .lt. 1) im=nx
      rhsphi(i,j)=rhsphi(i,j)+gamma*(((phi(ip,j)**2.d0-phi(ip,j))*normx(ip,j) - (phi(im,j)**2.d0-phi(im,j))*normx(im,j))*0.5*dxi + &
                                     ((phi(i,jp)**2.d0-phi(i,jp))*normy(i,jp) - (phi(i,jm)**2.d0-phi(i,jm))*normy(i,jm))*0.5*dyi)
    enddo
  enddo

  ! phase-field n+1 (Euler explicit)
  do i=1,nx
    do j=2,ny-1
      phi(i,j) = phi(i,j)  + dt*rhsphi(i,j);
    enddo
  enddo


  ! impose BC on the phase-field
  do i=1,nx
    phi(i,1) = phi(i,2)
    phi(i,ny) = phi(i,ny-1)
  enddo
  !$acc end kernels
  ! no flux at the walls
  !write(*,*) "maxphi", maxval(phi)
   write(*,*) "phi center", phi(32,32)


  !##########################################################
  ! STEP 1: Advance temperature field 
  !##########################################################
  ! Advection + diffusion (one loop, faster on GPU)
  !$acc kernels
  do i=1,nx
    do j=2,ny-1
      ip=i+1
      im=im-1
      jp=j+1
      jm=j-1
      if (ip .gt. nx) ip=1
      if (im .lt. 1) im=nx
      ! add advection contribution
      rhstemp(i,j)=-(u(ip,j)*0.5*(temp(ip,j)+temp(i,j)) - u(i,j)*0.5*(temp(i,j)+temp(im,j)))*dxi -(v(i,jp)*0.5*(temp(i,jp)+temp(i,j)) - v(i,j)*0.5*(temp(i,j)+temp(i,jm)))*dyi
      ! add diffusion contribution
      rhstemp(i,j)=rhstemp(i,j) + difftemp*((temp(ip,j)-2.d0*temp(i,j)+temp(im,j))*ddxi + (temp(i,jp) -2.d0*temp(i,j) +temp(i,jm))*ddyi)      
    enddo
  enddo

  ! New temperature field
  do i=1,nx
    do j=2,ny-1
      temp(i,j) = temp(i,j)  + dt*rhstemp(i,j);
    enddo
  enddo

  ! impose BC on the temperature field
  do i=1,nx
    temp(i,1) = 1.0d0
    temp(i,ny) = -1.0d0
  enddo
  !$acc end kernels

  !##########################################################
  ! 3B Start of Poisson solver, pressure in physical space obtained
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
  !open(unit=55,file='out.dat',form='unformatted',position='append',access='stream',status='new')
  !write(55) phi(:,:)
  !close(55)
  !##########################################################
  !End STEP 3B of Poisson solver, pressure in physical space obtained
  !##########################################################


  !##########################################################
  !START 3C: Start correction step
  !##########################################################
  !##########################################################
  !END 3C: End correction step
  !##########################################################
  call cpu_time(timef)
  print '(" Time elapsed = ",f6.1," ms")',1000*(timef-times)


enddo

deallocate(x,y)
!deallocate(a,b,c,d,sol)
deallocate(kx,kx2)
deallocate(rhsp,p,rhspc)

end program main
