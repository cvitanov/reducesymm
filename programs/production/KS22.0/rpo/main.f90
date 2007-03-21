program shootME

use nrtype
use ifc_newt
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
real(dp) :: T,kappa, ti,tf, h, h2
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

220 Format(F21.16)
221 Format(<d>F30.18)
222 Format(<d/2+1>F30.18)
223 Format(<4>F30.18)

allocate(v(d),vdum(d),bc(d))
allocate(a(d/2+1),adum(d/2+1),ai(d/2+1),af(d/2+1))
allocate(lin(d/2+1),f0(d/2+1),f1(d/2+1),f2(d/2+1),f3(d/2+1),e(d/2+1),e2(d/2+1))
allocate(f0dum(d/2+1),f1dum(d/2+1),f2dum(d/2+1),f3dum(d/2+1),edum(d/2+1),e2dum(d/2+1))
allocate(wR(d),wI(d))

open(19,file=trim(wd)//'/rpoGuess.dat')
 
	read(19,*) v(1:d)
 
close(19)

open(20,file=trim(wd)//'/periodsGuess.dat')
 
	read(20,*) T
	read(20,*) kappa
 
close(20)

call dfftw_plan_dft_r2c_1d(plan,d,v,a,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
a=a/size(v)

bc(1:d/2)=real(a(2:size(a)))
bc(d/2+1:d)= aimag(a(2:size(a)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call mnewtRPO(Ntrial,bc,tolbc,tolf,T,kappa,ksFJ_rpo)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print *,"T",T,"kappa",kappa

open(27,file=trim(wd)//'/periods.dat')
	write (27,220) T
	write (27,220) kappa
close(27)

if ( newton_condition_met .ne. 0) then
	if ( newton_condition_met .eq. 2) then
		ti=0.0_dp
		tf=T
		h=abs(tf-ti)/Nsteps
		h2=h/2.0_dp
		call SetLin_KS(lin)
		call etdrk4DiagPrefactors(lin,h,R,M,f0,f1,f2,f3,e,e2)
		call etdrk4DiagPrefactors(lin,h2,R,M,f0dum,f1dum,f2dum,f3dum,edum,e2dum)
		ti=0.0_dp
		ai=(0.0_dp,0.0_dp)
		ai(2:size(a))=bc(1:size(bc)/2)+ii*bc(size(bc)/2+1:size(bc))
		af=(0.0_dp,0.0_dp)
		Jac=UnitMatrix(d)
		call etdrk4DiagJDriverSh(ti,ai,Jac,h,Nsteps,tf,af,Jac,f0,f1,f2,f3,e,e2,f0dum,f1dum,f2dum,f3dum,edum,e2dum,Nplt,etdrk4diag,etdrk4DiagJhr,SetNlin_KS,SetANdiag_KS)
	end if
	wR=0.0_dp
	wI=0.0_dp
	Jac=matmul(Rr(kappa/L,d),Jac)
	call la_geesx(Jac,wR(1:d),wI(1:d),select=SelectLargeEig_r,sdim=sdim)
	print *,"eig", wR(1:sdim)+ii*wI(1:sdim)
	open(35,file=trim(wd)//'/Jdiag.dat')
	do i=1,sdim
		write(35,"(2F30.18)") wR(i)+ii*wI(i)
	enddo
	close(35)
else
	stop "Newton Condition didn't meet."
endif

ti=0.0_dp
ai=(0.0_dp,0.0_dp)
ai(2:size(a))=bc(1:size(bc)/2)+ii*bc(size(bc)/2+1:size(bc))
print *,sum(abs(ai))
af=(0.0_dp,0.0_dp)
tf=real(Nrep,dp)*T
h=abs(tf-ti)/(Nrep*Nsteps)
h2=h/2.0_dp

open (33,file=trim(wd)//'/timestep.dat')
write (33,220) h2
close(33)

call SetLin_KS(lin)
call etdrk4DiagPrefactors(lin,h2,R,M,f0,f1,f2,f3,e,e2)
call etdrk4DiagDriverS(ti,ai,2*Nrep*Nsteps,tf,af,f0,f1,f2,f3,e,e2,Nplt,SetNlin_KS)
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

a=(0.0_dp,0.0_dp)
a(2:size(a))=bc(1:size(bc)/2)+ii*bc(size(bc)/2+1:size(bc))

adum=a

call dfftw_plan_dft_c2r_1d(invplan,d,adum,v,FFTW_ESTIMATE)
call dfftw_execute(invplan)
call dfftw_destroy_plan(invplan)

open(28,file=trim(wd)//'/rpoUic.dat')
	write(28,221) v
close(28)


end program