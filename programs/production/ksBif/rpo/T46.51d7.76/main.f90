Program ksEquil

use nrtype
use ks
use ifc_newt
use ifc_integr
use ifc_util
use f95_lapack, only: LA_GEEV, LA_GEESX
USE LA_PRECISION, ONLY: WP => DP

implicit none

include "fftw3.f"

real(dp), dimension(:),allocatable :: v,vx,vxx,vxxx
complex(dpc), dimension(:),allocatable :: a,adum,af,ai,R_c,Raf
complex(dp), dimension(:), allocatable :: w
real(dp), dimension(:),allocatable :: bc,fvec
real(dp), dimension(:,:),allocatable :: fjac
complex(dp), dimension(:,:), allocatable :: fjacdum,vR
integer(i8b) :: invplan, plan ! needed by fftw3
real(dp), dimension(:), allocatable :: wR,wI
integer(i4b) :: k,i,j,sdim=8, Nf, Nsteps2
real(dp) :: dL, kappa, T,ti,tf,h,h2,fvmax=5d-2,fvmin=1d-6,sfvec,damp=5.0_dp
character*64 :: wd
integer(i4b) :: nargs
logical :: logicdum,fv

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
	read(21,*) dL
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
	read(21,*)
	read(21,*) Nf
close(21)

220 Format(F21.16)
221 Format(<d>F21.16)
222 Format(<d/2+1>F21.16)
! 

allocate(v(d),vx(d),vxx(d),vxxx(d),a(d/2+1),adum(d/2+1),af(d/2+1),Raf(d/2+1),bc(d),ai(d/2+1),w(d))
allocate( fvec(d),fjacdum(d,d), vR(d,d))
! allocate(v(d),vdum(d),bc(d))
!allocate(a(d/2+1),adum(d/2+1),ai(d/2+1),af(d/2+1))
allocate(lin(d/2+1),f0(d/2+1),f1(d/2+1),f2(d/2+1),f3(d/2+1),e(d/2+1),e2(d/2+1))
allocate(f0dum(d/2+1),f1dum(d/2+1),f2dum(d/2+1),f3dum(d/2+1),edum(d/2+1),e2dum(d/2+1))
allocate(wR(d),wI(d))

open(19,file=trim(wd)//'/ic.dat')
 
	read(19,*) v
 
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

open(35,file=trim(wd)//'/Jdiag.dat')
open(24,file=trim(wd)//'/equilU.dat')
open(26,file=trim(wd)//'/c_steady.dat')
open(27,file=trim(wd)//'/periods.dat')

print *, "T",T,"Nsteps",Nsteps,"Nplt",Nplt

do i=1,Nf

	print *, "Point #",i, "L=",L,"dL",dL
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call mnewtRPO(ntrial,bc,tolbc,tolf,T,kappa,ksFJ_rpo)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	print *,"T",T,"kappa",kappa
	

	write (27,220) T
	write (27,220) kappa

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
		write(35,"(<sdim>F20.10)") abs(wR(1:sdim)+ii*wI(1:sdim))
	else
		stop "Newton Condition didn't meet."
	endif

	a=(0,0)

	a(2:size(a))=bc(1:d/2)+ii*bc(d/2+1:d)

	call dfftw_plan_dft_c2r_1d(invplan,d,a,v,FFTW_ESTIMATE)
	call dfftw_execute(invplan)
	call dfftw_destroy_plan(invplan)

	write(24,221) v
	
	call FourierDif(v,vx,L,1)
	call FourierDif(v,vxx,L,2)
	call FourierDif(v,vxxx,L,3)
	
	write(26,"(2F21.16)") L, sum(v**2-vx-vxxx)/d
	print *,"E",sum(v**2-vx-vxxx)/d
	
	L=L+dL
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	fv=.false.
	do while(.not. fv)
		call SetLin_KS(lin)
		ti=0.0_dp
		tf=T
		h=abs(tf-ti)/Nsteps
		call etdrk4DiagPrefactors(lin,h,R,M,f0,f1,f2,f3,e,e2)
		ai=(0.0_dp,0.0_dp)
		ai(2:size(ai))=bc(1:size(bc)/2)+ii*bc(size(bc)/2+1:size(bc))
		af=(0.0_dp,0.0_dp)

		call etdrk4DiagDriverS(ti,ai,Nsteps,tf,af,f0,f1,f2,f3,e,e2,Nplt,SetNlin_KS)

		Raf=Rc(-kappa/L,size(af))*af
		fvec(1:d/2)=real(ai(2:size(ai))-Raf(2:size(ai)))
		fvec(d/2+1:d)=aimag(ai(2:size(ai))-Raf(2:size(ai)))

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 


		sfvec=sum(abs(fvec))
		print *,"sfvec",sfvec
		if ( sfvec .gt. fvmax ) then 
			L=L-dL
			dL=0.3_dp*dL
			L=L+dL
		else if ( sfvec .lt. fvmin ) then
			dL=2.0_dp*dL
			L=L+dL
		else
			fv=.true.
		endif
	enddo
	
end do

close(35)
close(26)
close(24)
close(27)

end program
