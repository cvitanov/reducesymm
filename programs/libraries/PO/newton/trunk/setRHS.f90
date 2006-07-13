subroutine setRHS(diff,RHS)

use nrtype
use nrutil, only:assert_eq

implicit none

real(dp),  dimension(:), intent(in) :: diff
real(dp),  dimension(:), intent(out) :: RHS

integer(i4b) :: ndum

ndum = assert_eq(size(diff)+1,size(RHS,1),'setRHS')

RHS=0.0_dp

RHS(1:size(diff)) = diff

end subroutine