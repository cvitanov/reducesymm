subroutine SetNlin_Roessler(a,N_a)

use nrtype
use parameters

implicit none

complex(dpc), dimension(:), intent(in) :: a
complex(dpc), dimension(:), intent(out) :: N_a

N_a(1) = -a(2)-a(3)
N_a(2) = a(1)
N_a(3) = beta+a(3)*a(1)

end subroutine
