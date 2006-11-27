program integrRPO

use nrtype
use ifc_integr
use f95_lapack, only: LA_GEESX
use la_precision, only: wp => dp
use ifc_util
use ks

implicit none

include "fftw3.f"

real(dp), dimension(:), allocatable :: v,vdum, bc
real(dp), dimension(:), allocatable :: wR,wI
complex(dpc), dimension(:),allocatable :: a,adum
complex(dpc), dimension(:), allocatable :: ai,af
integer(i8b) :: invplan, plan ! needed by fftw3
integer(i4b) :: k,i, sdim, Nrep=3
real(dp) :: T,kappa, ti=0.0_dp,tf=200.0_dp, h, h2
real(dp) :: damp=13.0_dp
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
	read(21,*) 
	read(21,*) Nsteps
	read(21,*) 
	read(21,*) Nplt
	read(21,*)
	read(21,*) M
	read(21,*)
	read(21,*) R
close(21)

220 Format(F30.18)
221 Format(<d>F30.18)
222 Format(<d/2+1>F30.18)
223 Format(<4>F30.18)

allocate(v(d),vdum(d),bc(d))
allocate(a(d/2+1),adum(d/2+1),ai(d/2+1),af(d/2+1))
allocate(lin(d/2+1),f0(d/2+1),f1(d/2+1),f2(d/2+1),f3(d/2+1),e(d/2+1),e2(d/2+1))
allocate(f0dum(d/2+1),f1dum(d/2+1),f2dum(d/2+1),f3dum(d/2+1),edum(d/2+1),e2dum(d/2+1))
allocate(wR(d),wI(d))

open(19,file=trim(wd)//'/rpoUic.dat')
 
	read(19,*) v(1:d)
 
close(19)

call dfftw_plan_dft_r2c_1d(plan,d,v,a,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
a=a/size(v)

bc(1:d/2)=real(a(2:size(a)))
bc(d/2+1:d)= aimag(a(2:size(a)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
bc=0.0_dp

open(21,file=trim(wd)//'tw.dat')
	do i=1,15
	read(21,*) bc(i)
	read(21,*) bc(i+1)
	end do
close(21)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



ti=0.0_dp
ai=(0.0_dp,0.0_dp)
ai(2:size(a))=bc(1:size(bc)/2)+ii*bc(size(bc)/2+1:size(bc))
print *,sum(abs(ai))
af=(0.0_dp,0.0_dp)

open (33,file=trim(wd)//'/timestep.dat')
read(33,220) h
close(33)

Nsteps=int(abs(tf-ti)/h,i4b)

tf=Nsteps*h

print *,"int",h,Nsteps,h*Nsteps,sum(abs(ai))

call SetLin_KS(lin)
call etdrk4DiagPrefactors(lin,h,R,M,f0,f1,f2,f3,e,e2)
call etdrk4DiagDriverS(ti,ai,Nsteps,tf,af,f0,f1,f2,f3,e,e2,Nplt,SetNlin_KS)
open (29,file=trim(wd)//'/rpoU.dat')
do i=1,size(aSt,1)
	adum=aSt(i,:)
!	print *,sum(abs(adum))
	call dfftw_plan_dft_c2r_1d(invplan,d,adum,v,FFTW_ESTIMATE)
	call dfftw_execute(invplan)
	call dfftw_destroy_plan(invplan)
	write(29,221) v		
end do
close(29)


end program