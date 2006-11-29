Program ksEquil

use nrtype
use ks
use ifc_newt
use ifc_util

implicit none

include "fftw3.f"

real(dp), dimension(:),allocatable :: v,vx,vxx
integer(i8b) :: invplan, plan ! needed by fftw3
integer(i4b) :: k,i,j
character*64 :: wd
integer(i4b) :: nargs

nargs=iargc()

if (nargs .ne. 1) then
	print *,"Program must be called with exactly one argument indicating the data directory."
	call exit
end if

call getarg(1,wd)

open(21,file=trim(wd)//'/parameters.dat')
	read(21,*)
	read(21,*) d
	read(21,*)
	read(21,*) L 
close(21)

220 Format(F30.18)
221 Format(<d>F21.16)
222 Format(<d/2+1>F21.16)
! 

allocate(v(d),vx(d),vxx(d))

open(19,file=trim(wd)//'/equilU.dat')
 
	read(19,*) v
 
close(19)

call FourierDif(v,vx,L,1)
call FourierDif(v,vxx,L,2)
 
open(24,file=trim(wd)//'/equilS.dat')

	write(24,221) v
	write(24,221) vx
	write(24,221) vxx

close(24)

end program