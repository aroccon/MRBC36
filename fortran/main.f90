program main
implicit none
integer, parameter :: nx=128, ny=128
double precision, parameter :: pi=141592653589793d0d0
double precision :: dx, dy, acoeff
integer :: i, j, k, n, m
double precision, allocatable :: x(:), y(:)
double precision, allocatable :: rhsp(:,:), y(:,:), p(:,:)

allocate(x(nx),y(ny))
allocate(rhsp(nx,ny),pext(nx,ny))

! Grid spacing: dx divided by nx to have perfect periodicity
dx = lx/nx 
dy = ly/(ny-1)
! Perturbation of test solution
n= 1
m= 2
acoeff= -1.d0/((n*n + (2*pi*m/ly)**2))

! Initialize rhsp term 
do j=1,ny
  do i=1,nx
    rhsp(i,j) = sin(n * x(i))*cos(2*pi*m * y(j)/ly)
    pext(i,j) = acoeff*rhsp(i,j)
  end do
end do


  ! FFT along x-direction
  do j=1,ny
    call fft_wrapper(f(:,j), f_hat(:,j), nx, 1)
  end do

  ! Solve for each wavenumber
  do k=1,nx
    ! Wavenumber for current mode
    if (k <= nx/2+1) then
      kx = 2.0d0*3.141592653589793d0*(k-1)/Lx
    else
      kx = 2.0d0*3.141592653589793d0*(k-1-nx)/Lx
    end if

    ! Tridiagonal coefficients (FD2 in y)
    do j=1,ny
      a(j) = 1.0d0/(dy**2)
      b(j) = -2.0d0/(dy**2) - kx**2
      c(j) = 1.0d0/(dy**2)
      rhs(j) = -real(f_hat(k,j))
    end do

    ! Apply Neumann BCs (no-flux)
    b(1) = -1.0d0/(dy**2) - kx**2  ! du/dy=0 at bottom
    c(1) = 1.0d0/(dy**2)
    a(ny) = 1.0d0/(dy**2)
    b(ny) = -1.0d0/(dy**2) - kx**2 ! du/dy=0 at top

    ! Gaussian elimination (no function calls)
    ! Forward elimination
    do j=2,ny
      b(j) = b(j) - a(j)/b(j-1)*c(j-1)
      rhs(j) = rhs(j) - a(j)/b(j-1)*rhs(j-1)
    end do

    ! Back substitution
    solution(ny) = rhs(ny)/b(ny)
    do j=ny-1,1,-1
      solution(j) = (rhs(j) - c(j)*solution(j+1))/b(j)
    end do

    u_hat(k,:) = cmplx(solution(:), 0.0d0)
  end do

  ! Inverse FFT
  do j=1,ny
    call fft_wrapper(u_hat(:,j), u(:,j), nx, -1)
  end do

  ! Output results
  open(unit=10, file='solution.dat')
  do j=1,ny
    do i=1,nx
      write(10,*) (i-1)*dx, (j-1)*dy, u(i,j)
    end do
  end do
  close(10)

contains

  ! Replace with actual FFT implementation
  subroutine fft_wrapper(in, out, n, dir)
    integer, intent(in) :: n, dir
    complex(8), intent(in) :: in(n)
    complex(8), intent(out) :: out(n)
    out = in  ! Placeholder - use FFTW or similar here
  end subroutine fft_wrapper

end program poisson_solver
