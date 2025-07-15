
program main

use cufft
use iso_c_binding
use openacc 
use velocity
use phase
use temperature

implicit none
integer, parameter :: nx=128, ny=64
double precision, parameter :: pi=3.141592653589793d0
double precision :: dx, dy, lx, ly, acoeff, q, l2norm, err, dyi, factor
integer :: i, j, k, n, m, t
integer :: ip,im,jp,jm
double precision, allocatable :: x(:), y(:), kx(:), kx2(:)
double complex :: a(ny), b(ny), c(ny), d(ny), sol(ny), meanp(ny)
integer :: planf, planb, status, ntmax, dump
double precision :: radius, eps, epsi, gamma, rho, mu, dxi, ddxi, ddyi, normod, dt, umax
double precision :: chempot, sigma
double precision :: pos, epsratio, times, timef, difftemp, h11, h12, h21, h22, rhoi, alpha, beta

!##########################################################
! declare parameters
! Define some basic quantities
ntmax=100
t=0
dump=10
radius=0.5d0
epsratio=1.0d0
gamma=0.0d0
sigma=0.0d0
radius=0.3d0
difftemp=7.0711e-04 ! sqrt(Ra) with Ra=2e6 
rho=1.d0
mu=7.0711e-04 ! sqrt(1/Ra) with Ra=2e6 
! Define domain size
lx = 2.d0
ly = 1.d0
! Grid spacing: dx divided by nx to have perfect periodicity, same for ny
dx = lx/nx 
dy = ly/ny
dxi=1.d0/dx
dyi=1.d0/dy
ddxi=dxi*dxi
ddyi=dyi*dyi
dt=0.0001
eps=max(dx,dy)
epsi=1.d0/eps
rhoi=1.d0/rho
alpha=1.d0
beta=0.0d0

!#define phiflag 0
!#define tempflag 0

!assign the code to one GPU
call acc_set_device_num(1,acc_device_nvidia)

!########################################################
! allocate variables (avoid allocate/deallocate at run time)
allocate(x(nx),y(ny))
allocate(rhsp(nx,ny),p(nx,ny),pext(nx,ny))
allocate(rhspc(nx/2+1,ny),pc(nx/2+1,ny))
allocate(kx(nx/2+1))
allocate(kx2(nx/2+1))

! velocity variables (defined on cell faces)
allocate(u(nx,ny),v(nx,ny+1))
allocate(rhsu(nx,ny),rhsv(nx,ny+1))
allocate(rhsu_o(nx,ny),rhsv_o(nx,ny+1))

! phase-field variables (defined on centers)
allocate(phi(nx,ny),rhsphi(nx,ny))
allocate(normx(nx,ny),normy(nx,ny))
allocate(fxst(nx,ny),fyst(nx,ny))

! temperature variables (defined on centers)
allocate(temp(nx,ny),rhstemp(nx,ny))


!########################################################

write(*,*) "NEMESI36"
write(*,*) "Code for 2D phase-field simulation in RB configuration"
write(*,*) "Grid.     :", nx, ny
write(*,*) "Lx and Ly :", lx, ly
write(*,*) "Dx and Dy  :", dx, dy
write(*,*) "Time step  :", dt
write(*,*) "eps        :", eps
write(*,*) "Density   :", rho
write(*,*) "Surf. ten :", sigma

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



!##########################################################
! Initialize fields
!##########################################################
write(*,*) "Initialize velocity, temperature and phase-field"
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
! temperature (internal + boundaries)
do j=2,ny-1
  do i=1,nx
    temp(i,j)=0.d0
  enddo
enddo
do i=1,nx
  temp(i,ny) = 0.0d0
  temp(i,1) =  1.0d0
enddo
! output fields
call writefield(t,1)
call writefield(t,2)
call writefield(t,3)
call writefield(t,4)
call writefield(t,5)
!##########################################################
! End fields init
!##########################################################




!##########################################################
! Start temporal loop
!##########################################################

write(*,*) "Start temporal loop"
do t=1,ntmax
  call cpu_time(times)
  write(*,*) "Time step",t,"of",ntmax
  !##########################################################
  ! START 1: Advance phase-field 
  !##########################################################
  ! Advection + diffusion term
  umax=0.0d0 ! replace then with real vel max
  gamma=1.0d0*umax
  do j=2,ny-1
    do i=1,nx
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

  ! compute sharpening (then move to ACDI)
  !$acc kernels
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
  !$acc end kernels

  ! phase-field n+1 (Euler explicit)
  !$acc kernels
  do i=1,nx
    do j=2,ny-1
      phi(i,j) = phi(i,j)  + dt*rhsphi(i,j);
    enddo
  enddo
  !$acc end kernels

  ! impose BC on the phase-field (no flux at the walls)
  !$acc kernels
  do i=1,nx
    phi(i,1) = phi(i,2)
    phi(i,ny) = phi(i,ny-1)
  enddo
  !$acc end kernels
  !write(*,*) "maxphi", maxval(phi)
  !write(*,*) "Phase field", phi(32,32)
  !##########################################################
  ! END 1: phase-field n+1 obtained
  !##########################################################



  !##########################################################
  ! STEP 2: Advance temperature field 
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
  ! Top wall hot and bottom wall cold
  do i=1,nx
    temp(i,1) =  1.0d0
    temp(i,ny) = 0.0d0
  enddo
  !$acc end kernels
  !##########################################################
  ! END 2: Temperature an n+1 obtained
  !##########################################################








  !##########################################################
  ! START 3A: Projection step for NS
  !##########################################################
  ! Advection + diffusion
  !$acc kernels
  do j=2,ny
    do i=1,nx
      ip=i+1
      im=im-1
      jp=j+1
      jm=j-1
      if (ip .gt. nx) ip=1
      if (im .lt. 1) im=nx
      ! compute the products (conservative form)
      h11 = (u(ip,j)+u(i,j))*(u(ip,j)+u(i,j))     - (u(i,j)+u(im,j))*(u(i,j)+u(im,j))
      h12 = (u(i,jp)+u(i,j))*(v(i,jp)+v(im,jp))   - (u(i,j)+u(i,jm))*(v(i,j)+v(im,j))
      h21 = (u(ip,j)+u(ip,jm))*(v(ip,j)+v(i,j))   - (u(i,j)+u(i,jm))*(v(i,j)+v(im,j))
      h22 = (v(i,jp)+v(i,j))*(v(i,jp)+v(i,j))     - (v(i,j)+v(i,jm))*(v(i,j)+v(i,jm))
      ! compute the derivative
      h11=h11*0.25d0*dxi
      h12=h12*0.25d0*dxi
      h21=h21*0.25d0*dxi
      h22=h22*0.25d0*dxi
      ! add advection to the rhs
      rhsu(i,j)=-(h11+h12)
      rhsv(i,j)=-(h21+h22)
      ! compute the diffusive terms
      h11 = mu*(u(ip,j)-2.d0*u(i,j)+u(im,j))*ddxi
      h12 = mu*(u(i,jp)-2.d0*u(i,j)+u(i,jm))*ddxi
      h21 = mu*(v(ip,j)-2.d0*v(i,j)+v(im,j))*ddxi
      h22 = mu*(v(i,jp)-2.d0*v(i,j)+v(i,jm))*ddxi
      rhsu(i,j)=rhsu(i,j)+(h11+h12)*rhoi
      rhsv(i,j)=rhsv(i,j)+(h21+h22)*rhoi
      ! add buoyancy term
      rhsv(i,j)=rhsv(i,j) + temp(i,j)+temp(i,jm)
    enddo
  enddo
  !$acc end kernels

  ! Compute surface tension forces 
  !$acc kernels
  do j=2,ny-1
    do i=1,nx
      ip=i+1
      im=i-1
      jp=j+1
      jm=j-1
      if (ip .gt. nx) ip=1
      if (im .lt. 1) im=nx
      chempot=phi(i,j)*(1.d0-phi(i,j))*(1.d0-2.d0*phi(i,j))*epsi-eps*(phi(ip,j)+phi(im,j)+phi(i,jp)+phi(i,jm)- 4.d0*phi(i,j))*ddxi
      fxst(i,j)=6.d0*sigma*chempot*0.5d0*(phi(ip,j)-phi(im,j))*dxi
      fyst(i,j)=6.d0*sigma*chempot*0.5d0*(phi(i,jp)-phi(i,jm))*dxi
    enddo
  enddo
  !$acc end kernels

  !$acc kernels
  ! Add surface tension forces to RHS (do not merge with above!)
  do j=2,ny-1
    do i=1,nx
      im=i-1
      jm=j-1
      if (im .lt. 1) im=nx
      !rhsu(i,j)=rhsu(i,j)+0.5d0*(fxst(im,j)+fxst(i,j))*rhoi
      !rhsv(i,j)=rhsv(i,j)+0.5d0*(fyst(i,jm)+fyst(i,j))*rhoi
    enddo
  enddo
  !$acc end kernels


  ! write surface tension forces (debug only)
 ! open(unit=55,file='out.dat',form='unformatted',position='append',access='stream',status='new')
 ! write(55) fyst(:,:)
 ! close(55)

  ! find u, v and w star (AB2), overwrite u,v and w
  !$acc kernels
  do j=2,ny
    do i=1,nx
      u(i,j) = u(i,j) + dt*(alpha*rhsu(i,j))!-beta*rhsu_o(i,j))
      v(i,j) = v(i,j) + dt*(alpha*rhsv(i,j))!-beta*rhsv_o(i,j))
      !rhsu_o(i,j)=rhsu(i,j)
      !rhsv_o(i,j)=rhsv(i,j)
    enddo
  enddo
  !$acc end kernels

  ! change after first loop AB2 coefficients
  alpha=1.5d0
  beta= 0.5d0

  !impose BCs on the flow field
  !$acc kernels
  do i=1,nx
    u(i,1)=0.d0
    u(i,ny)=0.0d0
    v(i,1)=0.0d0
    v(i,ny+1)=0.0d0
  enddo
  !$acc end kernels

  ! Compute rhs of Poisson equation div*ustar: divergence at cell center 
  !$acc kernels
  do j=1,ny-1
    do i=1,nx
      ip=i+1
      jp=j+1
      if (ip .gt. nx) ip=1
      rhsp(i,j) =             (rho*dxi/dt)*(u(ip,j)-u(i,j))
      rhsp(i,j) = rhsp(i,j) + (rho*dxi/dt)*(v(i,jp)-v(i,j))
    enddo
  enddo
  !$acc end kernels
  !##########################################################
  ! END 3A: Rhs (from projection computed)
  !##########################################################


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
    ! fix pressure on one point (otherwise is zero mean along y)
    if (i == 1) then
        b(1) = 1.0d0
        c(1) = 0.0d0
        d(1) = 0.0d0
    endif
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
  !write(55) p(:,:)
  !close(55)
  !##########################################################
  !End STEP 3B of Poisson solver, pressure in physical space obtained
  !##########################################################


  !##########################################################
  !START 3C: Start correction step
  !##########################################################
  !$acc kernels
  do j=1,ny-1
    do i=1,nx
      im=i-1
      jm=j-1
      if (im < 1) im=nx
      u(i,j)=u(i,j) - dt/rho*(p(i,j)-p(im,j))*dxi
      v(i,j)=v(i,j) - dt/rho*(p(i,j)-p(i,jm))*dxi
    enddo
  enddo

  ! re-impose BCs on the flow field
  do i=1,nx
    u(i,1)=0.d0
    u(i,ny)=0.0d0
    v(i,1)=0.0d0
    v(i,ny+1)=0.0d0
  enddo
  !$acc end kernels
  !##########################################################
  !END 3C: End correction step
  !##########################################################



  call cpu_time(timef)
  print '(" Time elapsed = ",f6.1," ms")',1000*(timef-times)

  !output fields
  if (mod(t,dump) .eq. 0) then
    write(*,*) "Saving output files"
    ! write velocity and pressure fiels (1-4)
	  call writefield(t,1)
	  call writefield(t,2)
	  call writefield(t,3)
    call writefield(t,4)
	  call writefield(t,5)
  endif

enddo

deallocate(x,y)
!deallocate(a,b,c,d,sol)
deallocate(kx,kx2)
deallocate(fxst,fyst,rhsu,rhsv,rhsu_o,rhsv_o)
deallocate(rhsp,p,rhspc)

end program main
