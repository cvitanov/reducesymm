subroutine setLHM(J,v,a,LHM)

use nrtype 
use nrutil, only: assert_eq
use ifc_util, only: UnitMatrix

implicit none

real(dp), dimension(:,:), intent(in) :: J
real(dp), dimension(:), intent(in) :: a, v
real(dp), dimension(:,:), intent(out) :: LHM 

integer(i4b) :: ndum,mdum

ndum = assert_eq(size(v),size(J,1),size(a),'setLHM 1')
mdum = assert_eq(ndum,size(LHM,1)-1,size(LHM,1)-1,'setLHM 2')

LHM(1:size(v),1:size(v)) = UnitMatrix(size(v)) - J
LHM(size(LHM,1),:) = a
LHM(:,size(LHM,1)) = v
LHM(size(LHM,1),size(LHM,1)) = 0.0_dp

end subroutine