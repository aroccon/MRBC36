subroutine readinput
use velocity
use phase
use param
implicit none

open(unit=55,file='input.inp',form='formatted',status='old')
!Time step parameters
read(55,*) restart
read(55,*) tstart
read(55,*) tfin
read(55,*) dump
!Flow parameters
read(55,*) dt
read(55,*) rho
read(55,*) mu
read(55,*) pr
read(55,*) lx
read(55,*) ly
! phase-field parameters
read(55,*) radius
read(55,*) sigma
read(55,*) epsr   



! compute pre-defined constant 
twopi=2.d0*pi 
difftemp=pr*mu
!difftemp=sqrt(pr/ra)! sqrt(Ra) with Ra=2e6 
!mu=sqrt(1.d0/ra) ! sqrt(1/Ra) with Ra=2e6 
dx = lx/(nx)
dy = ly/(ny-1)
dxi=1.d0/dx
ddxi=1.d0/dx/dx
dyi=1.d0/dy
ddyi=1.d0/dy/dy
rhoi=1.d0/rho
eps=epsr*min(dx,dy)
epsi=1.d0/eps
enum=1.e-16


write(*,*) "------------------------------------------------------"
write(*,*) "@@@@@@@@@@  @@@@@@@  @@@@@@@   @@@@@@@ @@@@@@    @@@@@"
write(*,*) "@@! @@! @@! @@!  @@@ @@!  @@@ !@@          @@! @@!@"    
write(*,*) "@!! !!@ @!@ @!@!!@!  @!@!@!@  !@!       @!!!:  @!@!@!@"
write(*,*) "!!:     !!: !!: :!!  !!:  !!! :!!          !!: !!:  !!!"
write(*,*) " :      :    :   : : :: : ::   :: :: : ::: ::   : : ::"
write(*,*) "------------------------------------------------------"
write(*,*) "Grid:    ", nx, 'x', ny
write(*,*) "Tfin     ", tfin
write(*,*) "Dump        ", dump
write(*,*) "Density      ", rho
write(*,*) "Viscosity      ", mu
write(*,*) "Prandtl        ", pr
write(*,*) "Radius          ", radius
write(*,*) "Sigma          ", sigma
write(*,*) "Eps             ", eps
write(*,*) 'Lx             ', lx
write(*,*) 'Ly             ', ly
write(*,*) 'Dx              ', dx
write(*,*) 'Dy              ', dy
end subroutine


