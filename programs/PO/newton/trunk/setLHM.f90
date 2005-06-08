subroutine setLHM(J,v,a,LHM)

use nrtype; use nrutil, only: assert_eq


real(dp), dimension(:,:), intent(in) :: J
real(dp), dimension(:), intent(in) :: a, v
real(dp), dimension(:,:), intent(out) :: LHM 

assert_eq(size(v),size(J,1),size(a),size(LHM,1)-1,,size(LHM,1)-1,'setLHM')

LHM(1:size(v),1:size(v)) = UnitMatrix(size(v),size(v)) - J
LHM(size(LHM,1),:) = a
LHM(:,size(LHM,1)) = v
LHM(size(LHM,1),size(LHM,1)) = 0.0_dp

end subroutine