module ks

use nrtype

real(dp), dimension(:,:), allocatable :: Jac
real(dp), dimension(:), allocatable :: lin,f0,f1,f2,f3,e,e2
real(dp), dimension(:), allocatable :: f0dum,f1dum,f2dum,f3dum,edum,e2dum
real(dp) :: R, L, Tw
integer(i4b) :: d, M, Nsteps, Nplt, Ntrial
real(dp) :: tolbc, tolf

interface ksFJ
	module procedure ksFJ_equil,ksFJ_rpo, ksFJ_req!(bc,fvec,fjac)
! 		use nrtype
! 		implicit none
! 		real(dp), DIMENSION(:), INTENT(IN) :: bc
! 		real(dp), DIMENSION(:), INTENT(OUT) :: fvec
! 		real(dp), DIMENSION(:,:), INTENT(OUT) :: fjac
! 	end subroutine
	!module subroutine ksFJ_rpo(bc,fvec,fjac,T,kappa)
! 		use nrtype
! 		implicit none
! 		real(dp), DIMENSION(:), INTENT(IN) :: bc
! 		real(dp), DIMENSION(:), INTENT(OUT) :: fvec
! 		real(dp), DIMENSION(:,:), INTENT(OUT) :: fjac
! 		real(dp), intent(in) :: T,kappa
! 	end subroutine
end interface ksFJ

! interface
! 	SUBROUTINE ksFJ_diff(bc,fvec,fjac,T,kappa)
! 		USE nrtype
! 		IMPLICIT NONE
! 		real(dp), DIMENSION(:), INTENT(IN) :: bc
! 		real(dp), DIMENSION(:), INTENT(OUT) :: fvec
! 		real(dp), DIMENSION(:,:), INTENT(OUT) :: fjac
! 		real(dp), intent(in) :: T,kappa
! 	end subroutine
! end interface


interface
	subroutine SetLin_KS(lnr)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(out) :: lnr
	end subroutine
end interface
interface
	subroutine SetNlin_KS(a,N_a)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in) :: a
		complex(dpc), dimension(:), intent(out) :: N_a
	end subroutine
end interface


INTERFACE
	SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE derivs	
END INTERFACE
interface
	function MatVar(x,y)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(size(y),size(y)) :: MatVar
	end function MatVar
end interface

interface Rc
	function Rc_M(w,dmn)
		use nrtype
		implicit none
		integer(i4b), intent(in) :: dmn
		complex(dpc), dimension(dmn) :: Rc_M
		real(dpc), intent(in) :: w
	end function Rc_M
	function Rc_a(w,a,dmn)
		use nrtype
		implicit none
		integer(i4b), intent(in) :: dmn
		complex(dpc), dimension(dmn) :: Rc_a
		complex(dpc), dimension(dmn), intent(in) :: a 
		real(dpc), intent(in) :: w
	end function Rc_a
end interface Rc

interface
	function Rr(w,dmn)
		use nrtype
		implicit none
		integer(i4b), intent(in) :: dmn
		real(dp), dimension(dmn,dmn) :: Rr
		real(dp), intent(in) :: w
	end function
end interface

interface
	subroutine SetANdiag_KS(a,Andiag)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in) :: a
		real(dp), dimension(:,:), intent(out) :: Andiag
	end subroutine
end interface

CONTAINS

SUBROUTINE ksFJ_equil(bc,fvec,fjac)
USE nrtype
use nrutil, only:assert_eq

IMPLICIT NONE

real(dp), DIMENSION(:), INTENT(IN) :: bc
real(dp), DIMENSION(:), INTENT(OUT) :: fvec
real(dp), DIMENSION(:,:), INTENT(OUT) :: fjac
!!!!
complex(dpc), dimension(size(bc)/2) :: a
complex(dpc), dimension(size(bc)/2+1) :: adum 
complex(dpc), dimension(size(bc)/2+1) :: N_adum
complex(dpc), dimension(size(bc)/2) :: fvec_c !, fvecA
real(dpc), dimension(d/2,d/2):: jcc, jbb, jbc, jcb
integer(i4b):: ndum,k,j
real(dp), dimension(size(bc)) :: v 
real(dp), dimension(size(bc)/2) :: q,linr
real(dp), dimension(size(fjac,1)):: wR,wI
complex(dpc), dimension(size(fjac,1)) :: W
complex(dpc), dimension(size(fjac,1),size(fjac,1)) :: VL, VR
integer(i4b) :: INFO

ndum=assert_eq(d,size(bc),size(fvec),'SetNlin1')
ndum=assert_eq(ndum,size(fjac,1),size(fjac,2),'SetNlin2')

! a does not include the a_0 coefficient
a=(0,0)
a=bc(1:d/2)+ ii*bc(d/2+1:d)
! 

adum=(0,0)
adum(2:size(adum))=a
call SetNlin_KS(adum,N_adum)

do k=1,d/2
	q(k)=k/L
	linr(k) = (1-(q(k))**2)*(q(k))**2
        fvec_c(k) = linr(k)*a(k) + N_adum(k+1)
end do


fvec(1:d/2)=real(fvec_c)
fvec(d/2+1:d)=aimag(fvec_c)

print *,"fvec", sum(abs(fvec))


jcc=0.0_dp
jbb=0.0_dp
jbc=0.0_dp
jcb=0.0_dp
fjac=0.0_dp

! Calculate Matrix of Variations(Jacobian)
!! calculate d\dot{c}/dc submatrix
do k=1,d/2
	jcc(k,k) = linr(k)
	do j=1,k-1
		jcc(k,j)=jcc(k,j)-2*q(k)*aimag(a(k-j))
	end do
	do j=k+1,d/2
		jcc(k,j)=jcc(k,j)+2*q(k)*aimag(a(j-k))
	end do
	do j=1,d/2-k
		jcc(k,j)=jcc(k,j)+2*q(k)*aimag(a(k+j))
	end do
end do
!! calculate d\dot{b}/db submatrix
do k=1,d/2
	jbb(k,k) = linr(k)
	do j=1,k-1
		jbb(k,j)=jbb(k,j)-2*q(k)*aimag(a(k-j))
	end do
	do j=k+1,d/2
		jbb(k,j)=jbb(k,j)+2*q(k)*aimag(a(j-k))
	end do
	do j=1,d/2-k
		jbb(k,j)=jbb(k,j)-2*q(k)*aimag(a(k+j))
	end do
end do
!! calculate d\dot{b}/dc submatrix
do k=1,d/2
	do j=1,k-1
		jbc(k,j)=-2*q(k)*real(a(k-j))
	end do
	do j=k+1,d/2
		jbc(k,j)=-2*q(k)*real(a(j-k))
	end do
	do j=1,d/2-k
		jbc(k,j)=jbc(k,j)+2*q(k)*real(a(k+j))
	end do
end do
!! calculate d\dot{c}/db submatrix
do k=1,d/2
	do j=1,k-1
		jcb(k,j)=2*q(k)*real(a(k-j))
	end do
	do j=k+1,d/2
		jcb(k,j)=jcb(k,j)+2*q(k)*real(a(j-k))
	end do
	do j=1,d/2-k
		jcb(k,j)=jcb(k,j)+2*q(k)*real(a(k+j))
	end do
end do

fjac(1:d/2,1:d/2)=jbb
fjac(1:d/2,d/2+1:d)=jbc
fjac(d/2+1:d,1:d/2)=jcb
fjac(d/2+1:d,d/2+1:d)=jcc

END SUBROUTINE ksFJ_equil

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ksFJ_rpo(bc,fvec,fjac,T,kappa)
USE nrtype
use nrutil, only:assert_eq
use ifc_integr
use ifc_util

IMPLICIT NONE

real(dp), DIMENSION(:), INTENT(IN) :: bc
real(dp), DIMENSION(:), INTENT(OUT) :: fvec
real(dp), DIMENSION(:,:), INTENT(OUT) :: fjac
real(dp), intent(in) :: T,kappa
!!!!
complex(dpc), dimension(size(bc)/2+1) :: ai,af,adum,Rai,Raf, R_cadum
complex(dpc), dimension(size(bc)/2+1) :: N_adum
complex(dpc), dimension(size(bc)/2+1) :: v_c,v_c_dum
real(dp), dimension(size(bc)) :: v,vi
integer(i4b):: k,i,ndum
real(dp) :: ti,tf,h,h2
complex(dpc), dimension(size(bc)/2+1) :: R_c, DR, DRRa
real(dp), dimension(size(bc),size(bc)) :: R_r,Ri_r

221 Format(<d>F21.16)
222 Format(<d/2+1>F21.16)
223 Format(<4>F21.16)

ndum=assert_eq(size(bc),size(fvec)-2,'ksFJ1')
ndum=assert_eq(d,size(fjac,1)-2,size(fjac,2)-2,'ksFJ2')

if ( .not. allocated(Jac)) allocate(Jac(d,d))
if ( .not. allocated(lin)) allocate(lin(d/2+1))

fvec=0.0_dp
ti=0.0_dp
tf=T

print *,"T",T
print *,"kappa",kappa

ai=(0.0_dp,0.0_dp)
ai(2:size(ai))= bc(1:d/2)+ ii*bc(d/2+1:d)
af=(0.0_dp,0.0_dp)

print *,"ai", sum(abs(ai))

Jac=UnitMatrix(d)

h=abs(tf-ti)/Nsteps

h2=h/2.0_dp

call SetLin_KS(lin)
call etdrk4DiagPrefactors(lin,h,R,M,f0,f1,f2,f3,e,e2)
call etdrk4DiagPrefactors(lin,h2,R,M,f0dum,f1dum,f2dum,f3dum,edum,e2dum)


call etdrk4DiagJDriverSh(ti,ai,Jac,h,Nsteps,tf,af,Jac,f0,f1,f2,f3,e,e2,f0dum,f1dum,f2dum,f3dum,edum,e2dum,Nplt,etdrk4diag,etdrk4DiagJhr,SetNlin_KS,SetANdiag_KS)

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

fjac(1:d,1:d)=UnitMatrix(d)-matmul(R_r,Jac)
fjac(1:d,d+1)=-v
fjac(1:d/2,d+2)=real(DRRa(2:d/2+1))
fjac(d/2+1:d,d+2)=aimag(DRRa(2:d/2+1))
fjac(d+2,1:d/2)=real(DR(2:d/2+1)*ai(2:d/2+1))
fjac(d+2,d/2+1:d)=aimag(DR(2:d/2+1)*ai(2:d/2+1))
fjac(d+1,1:d)= vi

print *,"fvec", sum(abs(fvec))

END SUBROUTINE ksFJ_rpo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE ksFJ_req(bc,fvec,fjac,kappa)
USE nrtype
use nrutil, only:assert_eq
use ifc_integr
use ifc_util

IMPLICIT NONE

real(dp), DIMENSION(:), INTENT(IN) :: bc
real(dp), DIMENSION(:), INTENT(OUT) :: fvec
real(dp), DIMENSION(:,:), INTENT(OUT) :: fjac
real(dp), intent(in) :: kappa
!!!!
complex(dpc), dimension(size(bc)/2+1) :: ai,af,adum,Rai,Raf, R_cadum
complex(dpc), dimension(size(bc)/2+1) :: N_adum
complex(dpc), dimension(size(bc)/2+1) :: v_c,v_c_dum
real(dp), dimension(size(bc)) :: v,vi
integer(i4b):: k,i,ndum
real(dp) :: ti,tf,h,h2
complex(dpc), dimension(size(bc)/2+1) :: R_c, DR, DRRa
real(dp), dimension(size(bc),size(bc)) :: R_r,Ri_r

221 Format(<d>F21.16)
222 Format(<d/2+1>F21.16)
223 Format(<4>F21.16)

ndum=assert_eq(size(bc),size(fvec)-1,'ksFJ1')
ndum=assert_eq(d,size(fjac,1)-1,size(fjac,2)-1,'ksFJ2')

if ( .not. allocated(Jac)) allocate(Jac(d,d))
if ( .not. allocated(lin)) allocate(lin(d/2+1))

fvec=0.0_dp
ti=0.0_dp
tf=Tw

print *,"T",Tw
print *,"kappa",kappa

ai=(0.0_dp,0.0_dp)
ai(2:size(ai))= bc(1:d/2)+ ii*bc(d/2+1:d)
af=(0.0_dp,0.0_dp)

print *,"ai", sum(abs(ai))

Jac=UnitMatrix(d)

h=abs(tf-ti)/Nsteps

h2=h/2.0_dp

call SetLin_KS(lin)
call etdrk4DiagPrefactors(lin,h,R,M,f0,f1,f2,f3,e,e2)
call etdrk4DiagPrefactors(lin,h2,R,M,f0dum,f1dum,f2dum,f3dum,edum,e2dum)


call etdrk4DiagJDriverSh(ti,ai,Jac,h,Nsteps,tf,af,Jac,f0,f1,f2,f3,e,e2,f0dum,f1dum,f2dum,f3dum,edum,e2dum,Nplt,etdrk4diag,etdrk4DiagJhr,SetNlin_KS,SetANdiag_KS)

R_c=Rc(-kappa/L,size(af))

R_r=Rr(kappa/L,d)

Raf=R_c*af

fvec(1:d/2)=real(ai(2:size(ai))-Raf(2:size(ai)))
fvec(d/2+1:d)=aimag(ai(2:size(ai))-Raf(2:size(ai)))
fvec(d+1)=0.0_dp

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

fjac(1:d,1:d)=UnitMatrix(d)-matmul(R_r,Jac)
fjac(1:d/2,d+1)=real(DRRa(2:d/2+1))
fjac(d/2+1:d,d+1)=aimag(DRRa(2:d/2+1))
fjac(d+1,1:d/2)=real(DR(2:d/2+1)*ai(2:d/2+1))
fjac(d+1,d/2+1:d)=aimag(DR(2:d/2+1)*ai(2:d/2+1))


print *,"fvec", sum(abs(fvec))

END SUBROUTINE ksFJ_req


end module