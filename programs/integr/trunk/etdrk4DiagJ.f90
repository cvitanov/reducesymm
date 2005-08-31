subroutine etdrk4Diag(a,h,aout,f0,f1,f2,f3,e,e2,SetNlin)

use nrtype
use nrutil, only: assert_eq

implicit none

complex(dpc), dimension(:), intent(in):: a ! Initial point
real(dp), intent(in) :: h ! Step size
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
! Advance solution by one step using etdrk4 algorithm (Diagonal case).
! Needs precomputed fuction f0-3,e,e2. SetNlin is calculated 
! at the various intermediate "positions".

integer(i4b) :: ndum
complex(dpc), dimension(size(a)) :: a1,a2,a3,Nlin_a,Nlin_a1,Nlin_a2,Nlin_a3 !tmp vars

ndum=assert_eq(size(f0), size(a), size(aout), 'etdrk4Diag-a' )

call setNlin(a,Nlin_a)
a1= e2*a + f0*Nlin_a 	! Take first step
J1= J + At(a)J
call setNlin(a1,Nlin_a1)
a2= e2*a + f0*Nlin_a1	! Take second step
J2= J + At(a1)J
call setNlin(a2,Nlin_a2)
a3= e2*a1 + f0*(2.0_dp*Nlin_a2-Nlin_a) ! Third step
call setNlin(a3,Nlin_a3)

aout = e*a + f1*Nlin_a + f2*(Nlin_a1+Nlin_a2) + f3*Nlin_a3 ! Sum them up approprietly 

end subroutine