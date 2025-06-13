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

do i=1,nx-1
  x(i+1)=x(i) + dx
enddo

write(*,*) "x", x

end program main
