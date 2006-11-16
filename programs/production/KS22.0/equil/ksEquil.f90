Program ksEquil

use nrtype
use ks
use ifc_newt
use ifc_util
use f95_lapack, only: LA_GEEV
USE LA_PRECISION, ONLY: WP => DP

implicit none

include "fftw3.f"

real(dp), dimension(:),allocatable :: v
complex(dpc), dimension(:),allocatable :: a,adum
complex(dp), dimension(:), allocatable :: w
real(dp), dimension(:),allocatable :: bc,fvec
real(dp), dimension(:,:),allocatable :: fjac
complex(dp), dimension(:,:), allocatable :: fjacdum,vR
real(dp), dimension(:),allocatable :: ar,ai
integer(i8b) :: invplan, plan ! needed by fftw3
integer(i4b) :: k,i,j,sdim=8
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

allocate(v(d),a(d/2+1),adum(d/2+1),bc(d),ar(d/2+1),ai(d/2+1),w(d),fvec(d),fjac(d,d),fjacdum(d,d), vR(d,d))

open(19,file=trim(wd)//'/equilGuess.dat')
 
	read(19,*) v
 
close(19)

 
call dfftw_plan_dft_r2c_1d(plan,d,v,a,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
a=a/size(v)

bc(1:d/2)=real(a(2:size(a)))
bc(d/2+1:d)= aimag(a(2:size(a)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call mnewt(ntrial,bc,tolbc,tolf,ksFJ_equil)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call ksFJ(bc,fvec,fjac)

fjacdum=fjac
w=0.0_dp
vR=0.0_dp
call la_geev(fjacdum,w,vR=vR) ! ,select=SelectSmallEig_c,sdim=sdim)
call sort_pick(w,vR)
print *,"eig", w(d-sdim:d)
open(35,file=trim(wd)//'/Jdiag.dat')
do i=d,d-sdim,-1
	write(35,"(2F30.18)") w(i)
enddo
close(35)

open(36,file=trim(wd)//'/Jev.dat')
do j=d,d-sdim,-1
	do i=1,d
		write(36,"(2F30.18)") vR(i,j)
	end do
enddo
close(36)

open(23,file=trim(wd)//'/equilPS.dat')

write(23,"(2F30.18)") 0.0_dp+ii*0.0_dp
do k=1,d/2
	write(23,"(2F30.18)") bc(k)+ii*bc(d/2+k)
end do

close(23)

a=(0,0)

a(2:size(a))=bc(1:d/2)+ii*bc(d/2+1:d)


call dfftw_plan_dft_c2r_1d(invplan,d,a,v,FFTW_ESTIMATE)
call dfftw_execute(invplan)
call dfftw_destroy_plan(invplan)

open(24,file=trim(wd)//'/equilU.dat')

write(24,221) v

close(24)

end program