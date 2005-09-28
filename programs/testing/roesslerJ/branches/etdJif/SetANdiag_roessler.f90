subroutine SetAndiag_roessler(a,ANdiag)

use nrtype
use parameters
implicit none

complex(dpc), dimension(:), intent(in) :: a
complex(dpc), dimension(:,:), intent(out) :: ANdiag

ANdiag(1,1) = (0.0,0.0)
ANdiag(1,2) = (-1.0,0.0)
ANdiag(1,3) = (-1.0,0.0)
ANdiag(2,1) = (1.0,0.0)
ANdiag(2,2) = (0.0,0.0)
ANdiag(2,3) = (0.0,0.0)
ANdiag(3,1) = a(3)
ANdiag(3,2) = (0.0,0.0)
ANdiag(3,3) = a(1)







end subroutine
