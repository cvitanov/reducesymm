SUBROUTINE ksFJ(bc,fvec,fjac,T,kappa)
USE nrtype
use nrutil, only:assert_eq
use ifc_integr
use ifc_util
use f95_lapack, only: LA_GEES
USE LA_PRECISION, ONLY: WP => DP
use ifc_rpo_ks

IMPLICIT NONE

include "fftw3.f"

real(dp), DIMENSION(:), INTENT(IN) :: bc
real(dp), DIMENSION(:), INTENT(OUT) :: fvec
real(dp), DIMENSION(:,:), INTENT(OUT) :: fjac
real(dp), intent(in) :: T,kappa
!real(dp), dimension(:), intent(in) :: q
!!!!
complex(dpc), dimension(size(bc)/2+1) :: ai,af,adum,Rai,Raf, R_cadum
complex(dpc), dimension(size(bc)/2+1) :: N_adum
complex(dpc), dimension(size(bc)/2+1) :: v_c,v_c_dum
real(dp), dimension(size(bc)) :: v,vi
! ! ! real(dp), dimension(size(bc)/2+1) :: Lin,f0,f1,f2,f3,e,e2
! ! ! real(dp), dimension(size(bc)/2+1) :: f0dum,f1dum,f2dum,f3dum,edum,e2dum
real(dp), dimension(size(bc),size(bc)):: Ji,Jf,Jdum
integer(i4b):: d,k,i
real(dp) :: ti,tf,h,h2
real(dp), dimension(size(fjac,1)):: wR,wI
integer(i8b) :: invplan, plan ! needed by fftw3
complex(dpc), dimension(size(bc)/2+1) :: R_c, DR, DRRa
real(dp), dimension(size(bc),size(bc)) :: R_r,Ri_r

221 Format(<d>F21.16)
222 Format(<d/2+1>F21.16)
223 Format(<4>F21.16)

d=assert_eq(size(bc),size(fvec)-2,'ksFJ1')
d=assert_eq(d,size(fjac,1)-2,size(fjac,2)-2,'ksFJ2')

fvec=0.0_dp
ti=0.0_dp
tf=T

print *,"T",T
print *,"kappa",kappa

ai=(0.0_dp,0.0_dp)
ai(2:size(ai))= bc(1:d/2)+ ii*bc(d/2+1:d)
af=(0.0_dp,0.0_dp)

Ji=UnitMatrix(d)

call SetLin_KS(lin)

h=abs(tf-ti)/Nsteps

h2=h/2.0_dp

call etdrk4DiagPrefactors(lin,h,R,M,f0,f1,f2,f3,e,e2) !! Move this to a module to reduce computational effort.

call etdrk4DiagPrefactors(lin,h2,R,M,f0dum,f1dum,f2dum,f3dum,edum,e2dum)

call etdrk4DiagJDriverSh(ti,ai,Ji,h,Nsteps,tf,af,Jf,f0,f1,f2,f3,e,e2,f0dum,f1dum,f2dum,f3dum,edum,e2dum,Nplt,etdrk4diag,etdrk4DiagJhr,SetNlin_KS,SetANdiag_KS)

open(17,file='a.dat')
open(21,file='v.dat')
open(27,file='aI.dat')
open(28,file='aR.dat')
	do k=1,size(aSt,1)	
		write(17,223) real(aSt(k,3)),aimag(aSt(k,3)),real(aSt(k,5)),tSt(k)
		write(27,222) aimag(aSt(k,:))
		write(28,222) real(aSt(k,:))
		adum=aSt(k,:)
		call dfftw_plan_dft_c2r_1d(invplan,d,adum,v,FFTW_ESTIMATE)
		call dfftw_execute(invplan)
		call dfftw_destroy_plan(invplan)
		write(21,221) v
	end do
!	write(17,'(4F15.10)') aimag(af(1))*L,aimag(af(2))*L,aimag(af(5))*L,tf/(L**2)	
close(17)
close(21)
close(27)
close(28)

! open(35,file='J.dat')
! 	do i=1,d
! 		write(35,221) Jf(i,:) 
! 	end do
! close(35)

!print *,"ai-af",ai-af

R_c=Rc(-kappa/L,size(af))

R_r=Rr(kappa/L,d)

Raf=R_c*af

fvec(1:d/2)=real(ai(2:size(ai))-Raf(2:size(ai)))
fvec(d/2+1:d)=aimag(ai(2:size(ai))-Raf(2:size(ai)))
fvec(d+1)=0.0_dp
fvec(d+2)=0.0_dp

!print *,"af",sum(abs(af)),sum(abs(Raf))
!print *,"fvec all", fvec

!!! Compute v at af for the LHS

call SetNlin_KS(af,N_adum)

do k=1,d/2+1
        v_c(k) = lin(k)*af(k) + N_adum(k)
end do

v_c=R_c*v_c

v(1:d/2)=real(v_c(2:d/2+1))
v(d/2+1:d)=aimag(v_c(2:d/2+1))

! !!! Compute v at R_c*af for comparison
! 
! call SetNlin_KS(Raf,N_adum)
! 
! do k=1,d/2+1
!         v_c_dum(k) = lin(k)*Raf(k) + N_adum(k)
! end do
! 
! print *,"compare", sum(abs(v_c_dum-v_c))

!!! Compute vi at ai for the LHS

call SetNlin_KS(ai,N_adum)

do k=1,d/2+1
        v_c(k) = lin(k)*ai(k) + N_adum(k)
end do

vi(1:d/2)=real(v_c(2:d/2+1))
vi(d/2+1:d)=aimag(v_c(2:d/2+1))

do k=1,d/2+1
	DR(k)=ii*Real((k-1),dp)/L
end do

DRRa= DR*(Raf)

fjac=0.0_dp

fjac(1:d,1:d)=UnitMatrix(d)-matmul(R_r,Jf)
fjac(1:d,d+1)=-v
fjac(1:d/2,d+2)=real(DRRa(2:d/2+1))
fjac(d/2+1:d,d+2)=aimag(DRRa(2:d/2+1))
fjac(d+2,1:d/2)=real(DR(2:d/2+1)*ai(2:d/2+1))
fjac(d+2,d/2+1:d)=aimag(DR(2:d/2+1)*ai(2:d/2+1))
fjac(d+1,1:d)= vi

Jdum=matmul(R_r,Jf)
!print *,"Jdum",Jdum(1,:)
! Jdum=Jf


wR=0.0_dp
wI=0.0_dp
call la_gees(Jdum,wR(1:d),wI(1:d))
print *,"eig", wR(1:d)+ii*wI(1:d)
print *,"fvec", sum(abs(fvec))

END SUBROUTINE
