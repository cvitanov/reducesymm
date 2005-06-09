subroutine setRHS(diff,RHS)

use nrtype
use nrutil, only:assert_eq

implicit none

real(dp),  dimension(:), intent(in) :: diff
real(dp),  dimension(:), intent(out) :: RHS

integer(i4b) :: ndum

ndum = assert_eq(size(diff),size(RHS,1),'setRHS')

RHS(1:size(RHS)-1) = diff
RHS(size(RHS)) = 0.0_dp 

end subroutine