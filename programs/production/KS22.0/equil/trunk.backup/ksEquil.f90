Program ksEquil

use nrtype
use parameters
use ifc_newt
use ifc_util
use f95_lapack, only: LA_GEESX
USE LA_PRECISION, ONLY: WP => DP

implicit none

include "fftw3.f"

real(dp), dimension(:),allocatable :: v
complex(dpc), dimension(:),allocatable :: a,adum
real(dp), dimension(:), allocatable :: wR,wI
real(dp), dimension(:),allocatable :: bc,fvec
real(dp), dimension(:,:),allocatable :: fjac
real(dp), dimension(:),allocatable :: ar,ai
integer(i8b) :: invplan, plan ! needed by fftw3
integer(i4b) :: k,i,sdim
!!!!
interface
	subroutine ksFJ(bc,fvec,fjac)
		USE nrtype
		implicit none
		real(dp), DIMENSION(:), INTENT(IN) :: bc
		real(dp), DIMENSION(:), INTENT(OUT) :: fvec
		real(dp), DIMENSION(:,:), INTENT(OUT) :: fjac
	end subroutine
end interface
!!!!
character*64 :: wd
integer(i4b) :: nargs
logical :: logicdum

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
	read(21,*)
	read(21,*) tolbc
	read(21,*)
	read(21,*) tolf
	read(21,*)
	read(21,*) Ntrial
close(21)

220 Format(F30.18)
221 Format(<d>F21.16)
222 Format(<d/2+1>F21.16)
! 

allocate(v(d),a(d/2+1),adum(d/2+1),bc(d),ar(d/2+1),ai(d/2+1),wR(d),wI(d),fvec(d),fjac(d,d))

open(19,file=trim(wd)//'/equilGuess.dat')
 
	read(19,*) v
 
close(19)

 
call dfftw_plan_dft_r2c_1d(plan,d,v,a,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
a=a/size(v)

bc(1:d/2)=real(a(2:size(a)))
bc(d/2+1:d)= aimag(a(2:size(a)))

call mnewt(ntrial,bc,tolbc,tolf,ksFJ)

call ksFJ(bc,fvec,fjac)


wR=0.0_dp
wI=0.0_dp
call la_geesx(fjac,wR,wI,select=SelectSmallEig,sdim=sdim)
print *,"eig", wR(1:sdim)+ii*wI(1:sdim)
open(35,file=trim(wd)//'/Jeig.dat')
do i=1,sdim
	write(35,"(2F30.18)") wR(i)+ii*wI(i)
enddo
close(35)

open(23,file=trim(wd)//'/equilR.dat')

write(23,*) 0.0_dp
do k=1,d/2
	write(23,*) bc(k)
end do

close(23)

open(25,file=trim(wd)//'/equilI.dat')

write(25,*) 0.0_dp
do k=d/2+1,d
	write(25,*) bc(k)
end do

close(23)

a=(0,0)

a(2:size(a))=bc(1:d/2)+ii*bc(d/2+1:d)


call dfftw_plan_dft_c2r_1d(invplan,d,a,v,FFTW_ESTIMATE)
call dfftw_execute(invplan)
call dfftw_destroy_plan(invplan)

open(24,file=trim(wd)//'/Uequil.dat')

write(24,221) v

close(24)

end program