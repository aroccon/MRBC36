
program main

use cufft
use param
use iso_c_binding
use openacc 
use velocity
use phase
use temperature
use nvtx

implicit none
integer :: i, j, k, n, m, t
integer :: ip,im,jp,jm
double precision, allocatable :: x(:), y(:), kx(:), kx2(:)
double precision, allocatable :: a(:,:), b(:,:), c(:,:) !TDMA pressure (real)
double complex,   allocatable :: d(:,:), sol(:,:) !TDMA pressure (complex variables)
double precision, allocatable :: af(:,:), bf(:,:), cf(:,:), df(:,:), solf(:,:) !TDMA temperature
integer :: planf, planb, status, stage
double precision, parameter ::  alpha(3)     = (/ 8.d0/15.d0,   5.d0/12.d0,   3.d0/4.d0 /) !rk3 alpha coef
double precision, parameter ::  beta(3)      = (/ 0.d0,       -17.d0/60.d0,  -5.d0/12.d0/) ! rk3 beta coef
!double precision, parameter ::  alpha(3) = (/ 1.d0,         3.d0/4.d0,    1.d0/3.d0 /) !rk3 ssp coef
!double precision, parameter ::  beta(3)  = (/ 0.d0,         1.d0/4.d0,    2.d0/3.d0 /) !rk3 ssp coef

#define phiflag 1
#define tempflag 1
#define impdifftemp 0

call readinput

!assign the code to one GPU
call acc_set_device_num(0,acc_device_nvidia)

!########################################################
! allocate variables (avoid allocate/deallocate at run time)
allocate(x(nx),y(ny))
allocate(rhsp(nx,ny),p(nx,ny))
allocate(rhspc(nx/2+1,ny),pc(nx/2+1,ny))
allocate(kx(nx/2+1))
allocate(kx2(nx/2+1))

! velocity variables (defined on cell faces)
allocate(u(nx,ny),v(nx,ny+1))
allocate(rhsu(nx,ny),rhsv(nx,ny+1))
allocate(rhsu_o(nx,ny),rhsv_o(nx,ny+1))
allocate(div(nx,ny))

! tdma variables (separated becasue they have different size)
allocate(a(nx/2+1,0:ny+1),b(nx/2+1,0:ny+1),c(nx/2+1,0:ny+1),d(nx/2+1,0:ny+1),sol(nx/2+1,0:ny+1)) ! pressure
#if impdifftemp == 1
allocate(af(nx,0:ny+1),bf(nx,0:ny+1),cf(nx,0:ny+1),df(nx,0:ny+1),solf(nx,0:ny+1)) ! temperature
#endif

! phase-field variables (defined on centers)
! add 0:ny+1, i.e. ghost nodes also for phi?
allocate(phi(nx,ny),rhsphi(nx,ny),psidi(nx,ny))
allocate(normx(nx,ny),normy(nx,ny))
allocate(fxst(nx,ny),fyst(nx,ny))

! temperature variables (defined on centers)
allocate(temp(nx,ny),rhstemp(nx,ny),rhstemp_o(nx,ny))
!########################################################


! Compute axis (position of cell centers)
x(1)=dx/2
y(1)=0.d0
do i=1,nx-1
  x(i+1)=x(i) + dx
enddo
do j=1,ny-1
  y(j+1)=y(j) + dy
enddo

! --- Wavenumbers and squares ---
do i = 1,nx/2+1
  kx(i) = (i-1)*(2.d0*pi/lx)
  kx2(i) = kx(i)**2
end do

! create plan forward x
status=0
status=status + cufftCreate(planf)
status=status + cufftPlanMany(planf, 1, nx, [nx, ny], 1, nx, [nx/2+1,ny], 1, nx/2+1, CUFFT_D2Z, ny) 
if (status.ne.0) write(*,*) "Error in cuFFT plan FWD:", status
! Create plan backward x
status=0
status=status + cufftCreate(planb)
status=status + cufftPlanMany(planb, 1, nx, [nx/2+1,ny], 1, nx/2+1, [nx,ny], 1, nx, CUFFT_Z2D, ny) 
if (status.ne.0) write(*,*) "Error in cuFFT plan BWD:", status


!##########################################################
! Initialize fields
!##########################################################
write(*,*) "Initialize velocity, temperature and phase-field"
! u velocity
do i=1,nx
  do j=1,ny
    u(i,j)= 0.1d0*sin(1.3d0*pi*(x(i)-dx/2))*cos(pi*(y(j)+dy/2))
  enddo
enddo
! v velocity
do i=1,nx
  do j=1,ny-1
    v(i,j)=-0.1d0*cos(pi*(x(i)))*sin(pi*(y(j)-dy/2))
  enddo
enddo
! phase-field
do i=1,nx
  do j=1,ny
    pos=(x(i)-lx/2)**2d0 + (y(j)-ly/2)**2d0
    phi(i,j)=0.5d0*(1.d0-tanh((sqrt(pos)-radius)/(2.d0*eps)))
  enddo
enddo
! temperature (internal + boundaries)
do j=2,ny-1
  do i=1,nx
    call random_number(noise)
    temp(i,j) = 1.0 - y(j) + 0.001d0*(2.0d0*noise - 1.0d0)
  enddo
enddo
do i=1,nx
  temp(i,ny) = 0.0d0
  temp(i,1) =  1.0d0
enddo
! output fields
call writefield(tstart,1)
call writefield(tstart,2)
call writefield(tstart,3)
#if phiflag == 1
call writefield(tstart,4)
#endif
# if tempflag == 1
call writefield(tstart,5)
#endif
!##########################################################
! End fields init
!##########################################################




!##########################################################
! Start temporal loop
!##########################################################
tstart=tstart+1
!$acc data copyin(u,v,phi,temp) copy(p) create(rhsu,rhsv,rhsphi,rhstemp,psidi,normx,normy,fxst,fyst,a,b,c,d,sol)
write(*,*) "Start temporal loop"
do t=tstart,tfin
  call cpu_time(times)
  write(*,*) "Time step",t,"of",tfin
  #if phiflag == 1
  !##########################################################
  ! START 1: Advance phase-field 
  !##########################################################
  ! Advection + diffusion term
  gamma=1.0d0*max(umax,vmax)
  write(*,*) "gamma", gamma
  !$acc kernels
  do j=2,ny-1
    do i=1,nx
      ip=i+1
      im=i-1
      jp=j+1
      jm=j-1
      if (ip .gt. nx) ip=1
      if (im .lt. 1) im=nx
      ! Advection + diffusion
      rhsphi(i,j)=-(u(ip,j)*0.5*(phi(ip,j)+phi(i,j)) - u(i,j)*0.5*(phi(i,j)+phi(im,j)))*dxi -(v(i,jp)*0.5*(phi(i,jp)+phi(i,j)) - v(i,j)*0.5*(phi(i,j)+phi(i,jm)))*dyi
      rhsphi(i,j)=rhsphi(i,j) + gamma*eps*((phi(ip,j)-2.d0*phi(i,j)+phi(im,j))*ddxi + (phi(i,jp) -2.d0*phi(i,j) +phi(i,jm))*ddyi) 
      val = max(0.d0, min(phi(i,j), 1.d0)) ! avoid machine precision overshoots in phi that leads to problem with log
      psidi(i,j) = eps*log((val+enum)/(1.d0-val+enum))
    enddo
  enddo

  ! compute normals from psidi
  do j=2,ny-1
    do i=1,nx
      ip=i+1
      im=i-1
      jp=j+1
      jm=j-1
      if (ip .gt. nx) ip=1
      if (im .lt. 1) im=nx
      ! Compute normals
      normx(i,j) = 0.5d0*(psidi(ip,j)-psidi(im,j))*dxi
      normy(i,j) = 0.5d0*(psidi(i,jp)-psidi(i,jm))*dyi
      normod = 1.0d0/(sqrt(normx(i,j)**2d0 + normy(i,j)**2d0) + enum)
      normx(i,j) = normx(i,j)*normod
      normy(i,j) = normy(i,j)*normod 
    enddo
  enddo
  ! gamma computed from previous field (after NS)
  do i=1,nx
    do j=2,ny-1
      ip=i+1
      im=i-1
      jp=j+1
      jm=j-1
      if (ip .gt. nx) ip=1
      if (im .lt. 1) im=nx
      rhsphi(i,j)=rhsphi(i,j)-gamma*((0.25d0*(1.d0-(tanh(0.5d0*psidi(ip,j)*epsi))**2)*normx(ip,j) - 0.25d0*(1.d0-(tanh(0.5d0*psidi(im,j)*epsi))**2)*normx(im,j))*0.5d0*dxi + &
                                     (0.25d0*(1.d0-(tanh(0.5d0*psidi(i,jp)*epsi))**2)*normy(i,jp) - 0.25d0*(1.d0-(tanh(0.5d0*psidi(i,jm)*epsi))**2)*normy(i,jm))*0.5d0*dyi) 
    enddo
  enddo

  ! phase-field n+1 (Euler explicit)
  do j=2,ny-1
    do i=1,nx
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
  #endif
  !write(*,*) "maxphi", maxval(phi)
  !write(*,*) "Phase field", phi(32,32)
  !##########################################################
  ! END 1: phase-field n+1 obtained
  !##########################################################



  !##########################################################
  ! STEP 2: Advance temperature field 
  !##########################################################
  #if tempflag == 1
  ! Advection + diffusion (one loop, faster on GPU)
  tempn=temp !tempn is the temperature at time step n (the initial one)
  do stage=1,3
    !$acc kernels
    do j=2,ny-1
      do i=1,nx
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
        if (ip .gt. nx) ip=1
        if (im .lt. 1) im=nx
        ! add advection contribution
        rhstemp(i,j)=-(u(ip,j)*0.5*(temp(ip,j)+temp(i,j)) - u(i,j)*0.5*(temp(i,j)+temp(im,j)))*dxi -(v(i,jp)*0.5*(temp(i,jp)+temp(i,j)) - v(i,j)*0.5*(temp(i,j)+temp(i,jm)))*dyi
        ! all diffusive contributions explicit
        #if impdifftemp == 0 
        rhstemp(i,j)=rhstemp(i,j) + difftemp*((temp(ip,j)-2.d0*temp(i,j)+temp(im,j))*ddxi + (temp(i,jp) -2.d0*temp(i,j) +temp(i,jm))*ddyi)   
        #endif
        ! x diffusive contribution explicit and y diffusive implicit (CN)
        #if impdifftemp == 1
        rhstemp(i,j)=rhstemp(i,j) + difftemp*(temp(ip,j)-2.d0*temp(i,j)+temp(im,j))*ddxi 
        #endif
      enddo
    enddo
    ! New provisional temperature field
    do j=1,ny
      do i=1,nx
        temp(i,j) = temp(i,j)  + dt*(alpha(stage)*rhstemp(i,j)+ alpha(stage)*rhstemp_o(i,j))
        rhstemp_o(i,j)=rhstemp(i,j)
      enddo
    enddo
    ! BC during RK stages
    do i=1,nx
      temp(i,1) =  1.0d0
      temp(i,ny) = 0.0d0
    enddo
    !$acc end kernels
  enddo

  #if impdifftemp == 1
  ! add the y diffusive of temperature implicit
  lambda = 0.5d0*difftemp*dt*ddyi ! then move in readinput
  !$acc parallel loop gang vector present(af, bf, cf, df, solf)
  do i=1,nx
    ! Build TDMA arrays for interior points only
    do j=2,ny-1
      af(i,j) = -lambda
      bf(i,j) =  1.0d0 + 2.0d0*lambda
      cf(i,j) = -lambda
      df(i,j) =  temp(i,j) + lambda*(temp(i,j+1)-2.0d0*temp(i,j) + temp(i,j-1))
    enddo
    ! Boundary conditions
    af(i,1)=0.0d0
    bf(i,1)=1.0d0
    cf(i,1)=0.0d0
    df(i,1)=1.0d0  ! bottom hot
    af(i,ny)=0.0d0
    bf(i,ny)=1.0d0
    cf(i,ny)=0.0d0
    df(i,ny)=0.0d0 ! top cold

    ! TDMA SOLVER - Forward sweep
    !$acc loop seq
    do j = 2, ny
      factor = af(i,j)/bf(i,j-1)
      bf(i,j) = bf(i,j) - factor*cf(i,j-1)
      df(i,j) = df(i,j) - factor*df(i,j-1)
    enddo
    ! Back substitution
    solf(i,ny) = df(i,ny)/bf(i,ny)
    !$acc loop seq
    do j=ny-1, 1, -1
      solf(i,j) = (df(i,j) - cf(i,j)*solf(i,j+1))/bf(i,j)
    end do
    do j=1, ny
      temp(i,j)=solf(i,j)
    enddo
  enddo
  #endif
  ! compute bottom and top nusselt numbers
  nut=0.0d0
  nub=0.0d0
  !$acc parallel loop reduction(+:nut,nub)
  do i=1,nx
    nut=nut + (temp(i,1)-temp(i,2))*dyi
    nub=nub + (temp(i,ny-1)-temp(i,ny))*dyi
  enddo
  nut=nut/nx
  nub=nub/nx
  #endif
  !##########################################################
  ! END 2: Temperature at n+1 obtained
  !##########################################################





  !##########################################################
  ! START 3A: Projection step for NS
  !##########################################################
  ! Advection + diffusion
  do stage=1,3
    !$acc kernels
    do j=2,ny-1
      do i=1,nx
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
        if (ip .gt. nx) ip=1
        if (im .lt. 1) im=nx
        ! compute the advection term (conservative form)
        h11 = 0.25d0*((u(ip,j)+u(i,j))*(u(ip,j)+u(i,j))     - (u(i,j)+u(im,j))*(u(i,j)+u(im,j)))*dxi
        h12 = 0.25d0*((u(i,jp)+u(i,j))*(v(i,jp)+v(im,jp))   - (u(i,j)+u(i,jm))*(v(i,j)+v(im,j)))*dyi
        h21 = 0.25d0*((u(ip,j)+u(ip,jm))*(v(ip,j)+v(i,j))   - (u(i,j)+u(i,jm))*(v(i,j)+v(im,j)))*dxi
        h22 = 0.25d0*((v(i,jp)+v(i,j))*(v(i,jp)+v(i,j))     - (v(i,j)+v(i,jm))*(v(i,j)+v(i,jm)))*dyi
        ! add advection to the rhs
        rhsu(i,j)=-(h11+h12)
        rhsv(i,j)=-(h21+h22)
        ! compute the diffusive terms
        h11 = mu*(u(ip,j)-2.d0*u(i,j)+u(im,j))*ddxi
        h12 = mu*(u(i,jp)-2.d0*u(i,j)+u(i,jm))*ddyi
        h21 = mu*(v(ip,j)-2.d0*v(i,j)+v(im,j))*ddxi
        h22 = mu*(v(i,jp)-2.d0*v(i,j)+v(i,jm))*ddyi
        rhsu(i,j)=rhsu(i,j)+(h11+h12)*rhoi
        rhsv(i,j)=rhsv(i,j)+(h21+h22)*rhoi
        ! add buoyancy term
        #if tempflag == 1
        rhsv(i,j)=rhsv(i,j) + 0.5d0*(temp(i,j) + temp(i,jm))
        #endif
        ! channel pressure driven (along x)
        !rhsu(i,j)=rhsu(i,j) + 1.d0
      enddo
    enddo
    !$acc end kernels

    ! Compute surface tension forces (LCSF mehthod, see Assessment of an energy-based surface tension model for simulation of...)
    #if phiflag == 1
    !$acc kernels
    do j=2,ny-1
      do i=1,nx
        ip=i+1
        im=i-1
        jp=j+1
        jm=j-1
        if (ip .gt. nx) ip=1
        if (im .lt. 1) im=nx
        curv=0.5d0*(normx(ip,j)-normx(im,j))*dxi + 0.5d0*(normy(i,jp)-normy(i,jm))*dyi
        fxst(i,j)= -sigma*curv*0.5d0*(phi(ip,j)-phi(im,j))*dxi!*phi(i,j)*(1.d0-phi(i,j))
        fyst(i,j)= -sigma*curv*0.5d0*(phi(i,jp)-phi(i,jm))*dyi!*phi(i,j)*(1.d0-phi(i,j))
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
        rhsu(i,j)=rhsu(i,j) + 0.5d0*(fxst(im,j)+fxst(i,j))*rhoi
        rhsv(i,j)=rhsv(i,j) + 0.5d0*(fyst(i,jm)+fyst(i,j))*rhoi
      enddo
    enddo
    !$acc end kernels
    #endif

    ! find u, v and w star (AB2), overwrite u,v and w
    !$acc kernels
    do j=2,ny-1
      do i=1,nx
        u(i,j) = u(i,j) + dt*(alpha(stage)*rhsu(i,j) + beta(stage)*rhsu_o(i,j))
        v(i,j) = v(i,j) + dt*(alpha(stage)*rhsv(i,j) + beta(stage)*rhsv_o(i,j))
        rhsu_o(i,j)=rhsu(i,j)
        rhsv_o(i,j)=rhsv(i,j)
      enddo
    enddo
    !$acc end kernels

    !impose BCs on the flow field
    !$acc kernels
    do i=1,nx
      u(i,1)=0.d0
      u(i,ny)=0.0d0
      v(i,1)=0.0d0
      v(i,ny+1)=0.0d0
    enddo
    !$acc end kernels

  enddo !end RK3 stage

  ! Compute rhs of Poisson equation div*ustar: divergence at cell center 
  !$acc kernels
  do j=1,ny
    do i=1,nx
      ip=i+1
      jp=j+1
      if (ip .gt. nx) ip=1
      rhsp(i,j) =             (rho*dxi/dt)*(u(ip,j)-u(i,j))
      rhsp(i,j) = rhsp(i,j) + (rho*dyi/dt)*(v(i,jp)-v(i,j))
    enddo
  enddo
  !$acc end kernels
  !##########################################################
  ! END 3A: Rhs (from projection computed)
  !##########################################################
  !call nvtxEndRange

  !call nvtxStartRange("Poisson")
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

  !$acc parallel loop gang vector present(a,b,c,d,sol,pc,rhspc,kx2)
  do i=1,nx/2+1
    ! Set up tridiagonal system for each kx
    ! The system is: (A_j) * pc(kx,j-1) + (B_j) * pc(kx,j) + (C_j) * pc(kx,j+1) = rhs(kx,j)
    ! FD2: -pc(j-1) + 2*pc(j) - pc(j+1)  --> Laplacian in y
    ! Neumann BC: d/dy pc = 0 at j=1 and j=ny
    ! Fill diagonals and rhs for each y
    do j = 1, ny
      a(i,j) =  1.0d0*ddyi
      b(i,j) = -2.0d0*ddyi - kx2(i)
      c(i,j) =  1.0d0*ddyi
      d(i,j) =  rhspc(i,j)
    end do
    ! Neumann BC at j=0 (ghost and first interior)
    b(i,0) = -1.0d0*ddyi - kx2(i)
    c(i,0) =  1.0d0*ddyi
    a(i,0) =  0.0d0
    ! Neumann BC at j=ny (top)
    a(i,ny+1) =  1.0d0*ddyi
    b(i,ny+1) = -1.0d0*ddyi - kx2(i)
    c(i,ny+1) =  0.0d0
    ! Special handling for kx=0 (mean mode)
    ! fix pressure on one point (otherwise is zero mean along x)
    if (i == 1) then
        b(i,1) = 1.0d0
        c(i,1) = 0.0d0
        d(i,1) = 0.0d0
    endif
    ! Thomas algorithm (TDMA) for tridiagonal system 
    ! Forward sweep
    do j = 1, ny+1
      factor = a(i,j)/b(i,j-1)
      b(i,j) = b(i,j) - factor * c(i,j-1)
      d(i,j) = d(i,j) - factor * d(i,j-1)
    end do
    ! Back substitution
    sol(i,ny+1) = d(i,ny+1) / b(i,ny+1)
    do j = ny, 0, -1
      sol(i,j) = (d(i,j) - c(i,j) * sol(i,j+1)) / b(i,j)
    end do

    ! Store solution (copy only interior nodes)
    do j = 1, ny
      pc(i,j) = sol(i,j)
    end do
  end do
  !!$acc end kernels

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
  !##########################################################
  !End STEP 3B of Poisson solver, pressure in physical space obtained
  !##########################################################
  
  !##########################################################
  !START 3C: Start correction step
  !##########################################################
  !$acc kernels
  do i=1,nx
    ! correct u
    do j=1,ny
      im=i-1
      if (im < 1) im=nx
      u(i,j)=u(i,j) - dt/rho*(p(i,j)-p(im,j))*dxi
    enddo
    do j=2,ny
      jm=j-1
      v(i,j)=v(i,j) - dt/rho*(p(i,j)-p(i,jm))*dyi
    enddo
  enddo
  !$acc end kernels

  cflx=umax*dt*dxi
  cfly=vmax*dt*dyi
  write(*,*) "CFL number:", max(cflx,cfly)
  write(*,*) "Mean nusselt", nut, nub


  ! re-impose BCs on the flow field
  umax=0.d0
  vmax=0.d0
  !$acc parallel loop collapse(1) reduction(+:umax,vmax)
  do i=1,nx
    u(i,1)=0.d0
    u(i,ny)=0.0d0
    v(i,1)=0.0d0
    v(i,ny+1)=0.0d0
    do j=2,ny
      umax=max(umax,u(i,j))
      vmax=max(vmax,v(i,j))
    enddo
  enddo

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
    #if phiflag == 1
    call writefield(t,4)
    #endif
    #if tempflag == 1
	  call writefield(t,5)
    #endif
  endif

enddo
!$acc end data

!write pressure (debug only)
!open(unit=55,file='output/div.dat',form='unformatted',position='append',access='stream',status='new')
!write(55) div
!close(55)

deallocate(x,y)
!deallocate(a,b,c,d,sol)
deallocate(kx,kx2)
deallocate(fxst,fyst,rhsu,rhsv,rhsu_o,rhsv_o)
deallocate(rhsp,p,rhspc)

end program main
