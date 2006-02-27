subroutine etdrk4DiagJ_g(a,J,h,aout,Jout,f0,f1,f2,f3,e,e2,SetNlin,SetAndiag)

use nrtype
use nrutil, only: assert_eq
use ifc_util

implicit none

complex(dpc), dimension(:), intent(in):: a ! Initial point
complex(dpc), dimension(:,:), intent(in):: J
real(dp), intent(in) :: h ! Step size
complex(dpc), dimension(:,:), intent(out):: Jout
real(dp), dimension(:),intent(in) :: f0,f1,f2,f3,e,e2 ! Functions of the linear operator
complex(dpc), dimension(:), intent(out):: aout ! Final point

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
		complex(dpc), dimension(:,:), intent(out) :: ANdiag
	end subroutine
end interface
! Advance solution by one step using etdrk4 algorithm (Diagonal case).
! Needs precomputed fuction f0-3,e,e2. SetNlin is calculated 
! at the various intermediate "positions".i

integer(i4b) :: d
complex(dpc), dimension(size(a)) :: a1,a2,a3,Nlin_a,Nlin_a1,Nlin_a2,Nlin_a3 !tmp vars
complex(dpc), dimension(size(a),size(a)) :: ANdiag
complex(dpc), dimension(size(a),size(a)) :: J1,J2,J3,Nlin_J,Nlin_J1,Nlin_J2,Nlin_J3


d=assert_eq(size(f0), size(a), size(aout), 'etdrk4Diag-a' )
d=assert_eq(size(a), size(J,1), size(Jout,1), 'etdrk4Diag-a' )

!Update $a$ first

call setNlin(a,Nlin_a)
a1= e2*a + f0*Nlin_a 	! Take first step
call setNlin(a1,Nlin_a1)
a2= e2*a + f0*Nlin_a1	! Take second step
call setNlin(a2,Nlin_a2)
a3= e2*a1 + f0*(2.0_dp*Nlin_a2-Nlin_a) ! Third step
call setNlin(a3,Nlin_a3)

aout = e*a + f1*Nlin_a + f2*(Nlin_a1+Nlin_a2) + f3*Nlin_a3 ! Sum them up approprietly 

!Then update \tilde{J} with points from the above calculation

call setAndiag(a,ANdiag)
Nlin_J = MatMul(ANdiag,J)
J1 = DiagMul(e2,J,d) + DiagMul(f0,Nlin_J,d)
call setAndiag(a1,ANdiag)
Nlin_J1 = MatMul(ANdiag,J1)
J2 = DiagMul(e2,J,d) + DiagMul(f0,Nlin_J1,d)
call setAndiag(a2,ANdiag)
Nlin_J2 = MatMul(ANdiag,J2)
J3 = DiagMul(e2,J1,d) + DiagMul(f0,2.0_dp*Nlin_J2-Nlin_J,d)
call setAndiag(a3,ANdiag)
Nlin_J3 = MatMul(ANdiag,J3)

Jout = DiagMul(e,J,d) + DiagMul(f1,Nlin_J,d) + DiagMul(f2,Nlin_J1+Nlin_J2,d) + DiagMul(f3,Nlin_J3,d)

end subroutine etdrk4DiagJ_g

subroutine etdrk4DiagJ_c0(a,J,h,aout,Jout,f0,f1,f2,f3,e,e2,SetNlin,SetAndiag,c0)

use nrtype
use nrutil, only: assert_eq
use ifc_util

implicit none

complex(dpc), dimension(:), intent(in):: a ! Initial point
complex(dpc), dimension(:,:), intent(in):: J
real(dp), intent(in) :: h ! Step size
complex(dpc), dimension(:,:), intent(out):: Jout
real(dp), dimension(:),intent(in) :: f0,f1,f2,f3,e,e2 ! Functions of the linear operator
complex(dpc), dimension(:), intent(out):: aout ! Final point
integer(i4b), intent(in), optional:: c0 !

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
		complex(dpc), dimension(:,:), intent(out) :: ANdiag
	end subroutine
end interface
! Advance solution by one step using etdrk4 algorithm (Diagonal case).
! Needs precomputed fuction f0-3,e,e2. SetNlin is calculated 
! at the various intermediate "positions".i

integer(i4b) :: d
complex(dpc), dimension(size(a)) :: a1,a2,a3,Nlin_a,Nlin_a1,Nlin_a2,Nlin_a3 !tmp vars
complex(dpc), dimension(size(a)-1,size(a)-1) :: ANdiag
complex(dpc), dimension(size(a)-1,size(a)-1) :: J1,J2,J3,Nlin_J,Nlin_J1,Nlin_J2,Nlin_J3


d=assert_eq(size(f0), size(a), size(aout), 'etdrk4Diag-a' )
d=assert_eq(size(a), size(J,1)+1, size(Jout,1)+1, 'etdrk4Diag-a' )

!Update $a$ first

call setNlin(a,Nlin_a)
a1= e2*a + f0*Nlin_a 	! Take first step
call setNlin(a1,Nlin_a1)
a2= e2*a + f0*Nlin_a1	! Take second step
call setNlin(a2,Nlin_a2)
a3= e2*a1 + f0*(2.0_dp*Nlin_a2-Nlin_a) ! Third step
call setNlin(a3,Nlin_a3)

aout = e*a + f1*Nlin_a + f2*(Nlin_a1+Nlin_a2) + f3*Nlin_a3 ! Sum them up approprietly 

!Then update \tilde{J} with points from the above calculation

call setAndiag(a,ANdiag)
Nlin_J = MatMul(ANdiag,J)
J1 = DiagMul(e2(2:size(a)),J,d-1) + DiagMul(f0(2:size(a)),Nlin_J,d-1)
call setAndiag(a1,ANdiag)
Nlin_J1 = MatMul(ANdiag,J1)
J2 = DiagMul(e2(2:size(a)),J,d-1) + DiagMul(f0(2:size(a)),Nlin_J1,d-1)
call setAndiag(a2,ANdiag)
Nlin_J2 = MatMul(ANdiag,J2)
J3 = DiagMul(e2(2:size(a)),J1,d-1) + DiagMul(f0(2:size(a)),2.0_dp*Nlin_J2-Nlin_J,d-1)
call setAndiag(a3,ANdiag)
Nlin_J3 = MatMul(ANdiag,J3)

Jout = DiagMul(e(2:size(a)),J,d-1) + DiagMul(f1(2:size(a)),Nlin_J,d-1) + DiagMul(f2(2:size(a)),Nlin_J1+Nlin_J2,d-1) + DiagMul(f3(2:size(a)),Nlin_J3,d-1)

end subroutine etdrk4DiagJ_c0
