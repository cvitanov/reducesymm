subroutine setRHS(diff,RHS)

use nrtype
use nrutil, only:assert_eq

implicit none

complex(dpc),  dimension(:), intent(in) :: diff
complex(dpc),  dimension(:), intent(out) :: RHS
!!
integer(i4b) :: ndum

ndum = assert_eq(size(diff)+2,size(RHS),'setRHS')

RHS=0.0_dp

RHS(1:size(diff)) = diff


end subroutine