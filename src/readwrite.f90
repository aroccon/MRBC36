
!##########################################################################
!###########################################################################
subroutine writefield(t,fieldn)
! Output field, file is written in the src/output folder 

use param
use velocity
use phase
use temperature 

implicit none
integer :: t,fieldn
character(len=40) :: namefile

! fieldn=1 means u
! fieldn=2 means v
! fieldn=3 means p
! fieldn=4 means phi
! fieldn=5 means temp

if (fieldn .eq. 1) then
write(namefile,'(a,i8.8,a)') './output/u_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
write(55) u
close(55)
endif

if (fieldn .eq. 2) then
write(namefile,'(a,i8.8,a)') './output/v_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
write(55) v
close(55)
endif

if (fieldn .eq. 3) then
write(namefile,'(a,i8.8,a)') './output/p_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
write(55) p
close(55)
endif

if (fieldn .eq. 4) then
write(namefile,'(a,i8.8,a)') './output/phi_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
write(55) phi
close(55)
endif

if (fieldn .eq. 5) then
write(namefile,'(a,i8.8,a)') './output/t_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',position='append',access='stream',status='new')
write(55) temp
close(55)
endif
end subroutine








subroutine readfield(t,fieldn)
! Used in case of fresh start, file is read from the src/init folder 
use velocity
use phase
use temperature
implicit none
integer :: t,fieldn
character(len=40) :: namefile


! fieldn=1 means u
! fieldn=2 means v
! fieldn=3 means p
! fieldn=4 means phi
! fieldn=5 means temp

if (fieldn .eq. 1) then
write(namefile,'(a,i8.8,a)') './init/u.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) u
close(55)
endif

if (fieldn .eq. 2) then
write(namefile,'(a,i8.8,a)') './init/v.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) v
close(55)
endif


if (fieldn .eq. 3) then
write(namefile,'(a,i8.8,a)') './init/p.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) p
close(55)
endif

if (fieldn .eq. 4) then
write(namefile,'(a,i8.8,a)') './init/phi.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) phi
close(55)
endif

if (fieldn .eq. 5) then
write(namefile,'(a,i8.8,a)') './init/t.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) temp
close(55)
endif
end subroutine






subroutine readfield_restart(t,fieldn)
! Used in case of restart, file is read from the src/output folder (iteration tstart must be present!)
use velocity
use phase
use temperature
implicit none
integer :: t,fieldn
character(len=40) :: namefile

! fieldn=1 means u
! fieldn=2 means v
! fieldn=4 means phi
! fieldn=5 means temp

if (fieldn .eq. 1) then
write(*,*) "Read u"
write(namefile,'(a,i8.8,a)') './output/u_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) u
close(55)
endif

if (fieldn .eq. 2) then
write(*,*) "Read v"
write(namefile,'(a,i8.8,a)') './output/v_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) v
close(55)
endif


if (fieldn .eq. 3) then
write(namefile,'(a,i8.8,a)') './output/p_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) p
close(55)
endif

if (fieldn .eq. 4) then
write(namefile,'(a,i8.8,a)') './output/phi_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) phi
close(55)
endif

if (fieldn .eq. 5) then
write(namefile,'(a,i8.8,a)') './output/t_',t,'.dat'
open(unit=55,file=namefile,form='unformatted',access='stream',status='old',convert='little_endian')
read(55) temp
close(55)
endif

end subroutine















subroutine nucheck(t)
! log the nussetl number
use param
use velocity
use phase
use temperature
implicit none
integer :: t,fieldn
character(len=40) :: namefile

! first iteration, create file
if (t .eq. 0) then
    open(unit=10, file='output/nu.dat', status='new', action='write')
    write(10, '(I10, 2F12.6)') t, nut, nub
    close(10)
endif

! random iteration, append on file
if (t .ne. 0) then
    open(unit=10, file='output/nu.dat', status='old', position='append', action='write')
    write(10, '(I10, 2F12.6)') t, nut, nub
    close(10)
endif

end subroutine
