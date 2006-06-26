subroutine etdrk4DiagJhr(a,adum,af,J,h,Jout,f0,f1,f2,f3,e,e2,SetNlin,SetAndiag)

use nrtype
use nrutil, only: assert_eq
use ifc_util

implicit none

complex(dpc), dimension(:), intent(in):: a,adum,af ! Initial,middle and final points
real(dp), dimension(:,:), intent(in):: J
real(dp), intent(in) :: h ! Step size
real(dp), dimension(:,:), intent(out):: Jout
real(dp), dimension(:),intent(in) :: f0,f1,f2,f3,e,e2 ! Functions of the linear operator
interface
	subroutine SetNlin(a,N_a)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in) :: a
		complex(dpc), dimension(:), intent(out) :: N_a
	end subroutine
end interface
interface
	subroutine SetANdiag(a,ANdiag)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in) :: a
		real(dp), dimension(:,:), intent(out) :: ANdiag
	end subroutine
end interface
! Advance the Jacobian by one step using etdrk4 algorithm (Diagonal case) for 
! a system with a Matrix of Variations computed in a real basis. 
! Needs precomputed fuction f0-3,e,e2. 

integer(i4b) :: dmn,dmnA
complex(dpc), dimension(size(a)) :: a2,a3,Nlin_a,Nlin_a1,Nlin_a2,Nlin_a3 !tmp vars
real(dp), dimension(2*(size(a)-1),2*(size(a)-1)) :: ANdiag
real(dp), dimension(2*(size(a)-1),2*(size(a)-1)) :: J1,J2,J3,Nlin_J,Nlin_J1,Nlin_J2,Nlin_J3
real(dp), dimension(2*(size(a)-1)) :: f0A,f1A,f2A,f3A,eA,e2A ! Functions of the linear operator to be used in calculation of Jacobian
complex(dpc), dimension(size(a)-1):: aA,a2A,afA

dmn=assert_eq(size(f0), size(a), size(af), 'etdrk4DiagJ-a' )
dmnA=assert_eq(2*(size(a)-1), size(J,1), size(Jout,1), 'etdrk4DiagJ-b' )

f0A=0.0_dp
f1A=0.0_dp
f2A=0.0_dp
f3A=0.0_dp
eA=0.0_dp
e2A=0.0_dp

aA=a(2:dmn)
a2A=adum(2:dmn)
afA=af(2:dmn)
f0A(1:dmn-1)=f0(2:dmn)
f0A(dmn:dmnA)=f0(2:dmn)
f1A(1:dmn-1)=f1(2:dmn)
f1A(dmn:dmnA)=f1(2:dmn)
f2A(1:dmn-1)=f2(2:dmn)
f2A(dmn:dmnA)=f2(2:dmn)
f3A(1:dmn-1)=f3(2:dmn)
f3A(dmn:dmnA)=f3(2:dmn)
eA(1:dmn-1)=e(2:dmn)
eA(dmn:dmnA)=e(2:dmn)
e2A(1:dmn-1)=e2(2:dmn)
e2A(dmn:dmnA)=e2(2:dmn)

call setAndiag(aA,ANdiag)
Nlin_J = MatMul(ANdiag,J)
J1 = DiagMul(e2A,J,dmnA) + DiagMul(f0A,Nlin_J,dmnA)
call setAndiag(a2A,ANdiag)
Nlin_J1 = MatMul(ANdiag,J1)
J2 = DiagMul(e2A,J,dmnA) + DiagMul(f0A,Nlin_J1,dmnA)
call setAndiag(a2A,ANdiag)
Nlin_J2 = MatMul(ANdiag,J2)
J3 = DiagMul(e2A,J1,dmnA) + 2.0_dp*DiagMul(f0A,Nlin_J2,dmnA)-DiagMul(f0A,Nlin_J,dmnA)	!+ DiagMul(f0A,2.0_dp*Nlin_J2-Nlin_J,dmnA)
call setAndiag(afA,ANdiag)
Nlin_J3 = MatMul(ANdiag,J3)

Jout = DiagMul(eA,J,dmnA) + DiagMul(f1A,Nlin_J,dmnA) + DiagMul(f2A,Nlin_J1,dmnA) +DiagMul(f2A,Nlin_J2,dmnA)+ DiagMul(f3A,Nlin_J3,dmnA)

end subroutine etdrk4DiagJhr

