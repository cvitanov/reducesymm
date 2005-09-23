subroutine SetLin_Roessler(Lin)

use nrtype
use parameters

implicit none

real(dp), dimension(:), intent(out) :: Lin

Lin=(/ 0.0_dp, Real(alpha), -Real(gamma) /)


end subroutine