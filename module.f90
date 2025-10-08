module param
  integer, parameter :: nx=512, ny=256
  double precision, parameter :: pi=3.141592653589793d0
  double precision :: dx, dy, lx, ly, acoeff, q, l2norm, err, dyi, factor, twopi
  double precision :: radius, eps, epsi, gamma, rho, mu, dxi, ddxi, ddyi, normod, dt
  double precision :: umax=0.0d0, vmax=0.d0, val, lambda
  double precision :: chempot, curv, sigma, cflx, cfly, ra, pr, nut, nub, num, noise, enum
  double precision :: pos, epsr, times, timef, difftemp, h11, h12, h21, h22, rhoi, alphag
  double precision :: ttop, tbot
  integer :: tstart, tfin, restart, dump, icphi
end module param

module velocity
    double precision, allocatable :: u(:,:),  v(:,:), rhsu(:,:), rhsv(:,:), rhsu_o(:,:), rhsv_o(:,:)
    double precision, allocatable :: rhsp(:,:),  p(:,:), div(:,:)
    double complex, allocatable :: rhspc(:,:), pc(:,:), rhs(:)
end module velocity

module phase
    double precision, allocatable :: rhsphi(:,:), phi(:,:), psidi(:,:), normx(:,:), normy(:,:), fxst(:,:), fyst(:,:)
end module phase

module temperature
    double precision, allocatable :: temp(:,:), tempn(:,:), rhstemp(:,:), rhstemp_o(:,:)
end module temperature


!module nvtx
!use iso_c_binding
!implicit none
!integer,private :: col(7) = [ int(Z'0000ff00'), int(Z'000000ff'), int(Z'00ffff00'), int(Z'00ff00ff'), int(Z'0000ffff'), int(Z'00ff0000'), int(Z'00ffffff')]
!character,private,target :: tempName(256)

!type, bind(C):: nvtxEventAttributes
!  integer(C_INT16_T):: version=1
!  integer(C_INT16_T):: size=48 !
!  integer(C_INT):: category=0
!  integer(C_INT):: colorType=1 ! NVTX_COLOR_ARGB = 1
!  integer(C_INT):: color
!  integer(C_INT):: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
!  integer(C_INT):: reserved0
!  integer(C_INT64_T):: payload   ! union uint,int,double
!  integer(C_INT):: messageType=1  ! NVTX_MESSAGE_TYPE_ASCII     = 1 
!  type(C_PTR):: message  ! ascii char
!end type

!interface nvtxRangePush
!  ! push range with custom label and standard color
!  subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
!  use iso_c_binding
!  character(kind=C_CHAR) :: name(256)
!  end subroutine

!  ! push range with custom label and custom color
!  subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
!  use iso_c_binding
!  import:: nvtxEventAttributes
!  type(nvtxEventAttributes):: event
!  end subroutine
!end interface

!interface nvtxRangePop
!  subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
!  end subroutine
!end interface

!contains

!subroutine nvtxStartRange(name,id)
!  character(kind=c_char,len=*) :: name
!  integer, optional:: id
!  type(nvtxEventAttributes):: event
!  character(kind=c_char,len=256) :: trimmed_name
!  integer:: i

!  trimmed_name=trim(name)//c_null_char

!  ! move scalar trimmed_name into character array tempName
!  do i=1,LEN(trim(name)) + 1
!     tempName(i) = trimmed_name(i:i)
!  enddo


!  if ( .not. present(id)) then
!    call nvtxRangePush(tempName)
!  else
!!    event%color=col(mod(id,7)+1)
!    event%message=c_loc(tempName)
!    call nvtxRangePushEx(event)
!  end if
!end subroutine

!subroutine nvtxEndRange
!  call nvtxRangePop
!end subroutine

!end module nvtx