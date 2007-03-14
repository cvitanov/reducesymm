Program ksEquil

use nrtype
use ks
use ifc_newt
use ifc_util
use f95_lapack, only: LA_GEEV
USE LA_PRECISION, ONLY: WP => DP

implicit none

include "fftw3.f"

real(dp), dimension(:),allocatable :: v,vx,vxx,vxxx
complex(dpc), dimension(:),allocatable :: a,adum
complex(dp), dimension(:), allocatable :: w
real(dp), dimension(:),allocatable :: bc,fvec
real(dp), dimension(:,:),allocatable :: fjac
complex(dp), dimension(:,:), allocatable :: fjacdum,vR
real(dp), dimension(:),allocatable :: ar,ai
integer(i8b) :: invplan, plan ! needed by fftw3
integer(i4b) :: k,i,j, sdim=10
character*64 :: wd
integer(i4b) :: nargs, elim=1
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

220 Format(F21.16)
221 Format(<d>F21.16)
222 Format(<d/2+1>F21.16)



allocate(v(d),vx(d),vxx(d),vxxx(d),a(d/2+1),adum(d/2+1),bc(d),ar(d/2+1),ai(d/2+1),w(d))
allocate( fvec(d),fjac(d,d),fjacdum(d,d), vR(d,d))

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
call mnewt_elim(ntrial,elim,bc,tolbc,tolf,ksFJ_equil_elim)
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
open(37,file=trim(wd)//'/JevU.dat')
do j=d,d-sdim,-1
	do i=1,d
		write(36,"(2F30.18)") vR(i,j)
	end do
	call bc2u(real(vR(:,j)),v)
	write(37,221) v
	call bc2u(aimag(vR(:,j)),v)
	write(37,221) v
enddo
close(36)
close(37)


adum=(0.0,0.0)
adum(2:size(a))=vR(1:d/2,d)+ii*vR(d/2+1:d,d)

call dfftw_plan_dft_c2r_1d(invplan,d,adum,v,FFTW_ESTIMATE)
call dfftw_execute(invplan)
call dfftw_destroy_plan(invplan)

open(29,file=trim(wd)//'/ev1.dat')
write(29,221) v
close(29)

adum=(0.0,0.0)
adum(2:size(a))=vR(1:d/2,d-1)+ii*vR(d/2+1:d,d-1)

call dfftw_plan_dft_c2r_1d(invplan,d,adum,v,FFTW_ESTIMATE)
call dfftw_execute(invplan)
call dfftw_destroy_plan(invplan)

open(29,file=trim(wd)//'/ev2.dat')
write(29,221) v
close(29)

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

call FourierDif(v,vx,L,1)
call FourierDif(v,vxx,L,2)
call FourierDif(v,vxxx,L,3)

open(25,file=trim(wd)//'equilS.dat')

	write(25,221) v(1:d)
	write(25,221) vx(1:d)
	write(25,221) vxx(1:d)

close(25)

print *,"vxxx", sum(v(1:d)),sum(v(1:d)**2), sum(vxxx(1:d)),sum(vxx(1:d)),sum(vx(1:d))
print *,"c",v(1)**2-vx(1)-vxxx(1),v(1),vx(1),vxxx(1)


open(26,file=trim(wd)//'steady_c.dat')
	write(26,220) sum(v**2-vx-vxxx)/d
close(25)

end program
