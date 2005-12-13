subroutine setLHM(f,J,v,q,LHM)

use nrtype 
use nrutil, only: assert_eq
use ifc_util, only: UnitMatrix, DiagMul

implicit none

complex(dpc), dimension(:,:), intent(in) :: J
complex(dpc), dimension(:), intent(in) ::  f, v, q
complex(dpc), dimension(:,:), intent(out) :: LHM 
!!
integer(i4b) :: ndum

ndum = assert_eq(size(v),size(J,1),size(f),'setLHM 1')
ndum = assert_eq(ndum,size(q),size(LHM,1)-1,'setLHM 2')
ndum = assert_eq(ndum,size(LHM,2)-1,'setLHM 3')


LHM=0.0_dp

LHM(1:size(v),1:size(v)) = UnitMatrix(size(v)) - J 

LHM(size(LHM,1),1:size(q)) = q

LHM(1:size(v),size(LHM,2)) = v

LHM(size(LHM,1),size(LHM,2)) = (1.0,0.0)

end subroutine