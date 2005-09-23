subroutine etdrk4DiagJ(a,J,h,aout,Jout,f0,f1,f2,f3,e,e2,SetNlin,SetAndiag)

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

integer(i4b) :: ndum
complex(dpc), dimension(size(a)) :: a1,a2,a3,Nlin_a,Nlin_a1,Nlin_a2,Nlin_a3 !tmp vars
complex(dpc), dimension(size(a),size(a)) :: ANdiag
complex(dpc), dimension(size(a),size(a)) :: kJ1,kJ2,kJ3,kJ4


ndum=assert_eq(size(f0), size(a), size(aout), 'etdrk4Diag-a' )

!Update $a$ first

call setNlin(a,Nlin_a)
a1= e2*a + f0*Nlin_a 	! Take first step
call setNlin(a1,Nlin_a1)
a2= e2*a + f0*Nlin_a1	! Take second step
call setNlin(a2,Nlin_a2)
a3= e2*a1 + f0*(2.0_dp*Nlin_a2-Nlin_a) ! Third step
call setNlin(a3,Nlin_a3)

aout = e*a + f1*Nlin_a + f2*(Nlin_a1+Nlin_a2) + f3*Nlin_a3 ! Sum them up approprietly 

!Then update \tilde{J} using IFRK4 with points from the above calculation

call setAndiag(a,ANdiag)
kJ1 = h*MatMul(ANdiag,J)
call setAndiag(a1,ANdiag)
kJ2 = h*MatMul(ANdiag,DiagMul(e2,J+kJ1/2,size(e2)) )
kJ3 = h*MatMul(ANdiag,DiagMul(e2,J,size(e2))+kJ2/2 )
call setANdiag(aout,ANdiag)
kJ4 = h*MatMul(ANdiag,DiagMul(e,J,size(e))+DiagMul(e2,kJ3,size(e2)))

Jout = DiagMul(e,J,size(e)) + DiagMul(e,kJ1,size(e))/6 + DiagMul(e2,kJ2 + kJ3,size(e2))/3 + kJ4/6

end subroutine
