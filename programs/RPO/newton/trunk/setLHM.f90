subroutine setLHM(f,J,v,q,R,diagk,LHM)

use nrtype 
use nrutil, only: assert_eq
use ifc_util, only: UnitMatrix, DiagMul

implicit none

complex(dpc), dimension(:,:), intent(inout) :: J
complex(dpc), dimension(:), intent(in) ::  f, v, q, R
real(dp), dimension(:), intent(in) :: diagk
complex(dpc), dimension(:,:), intent(out) :: LHM 
!!
integer(i4b) :: ndum

ndum = assert_eq(size(v),size(J,1),size(f),'setLHM 1')
ndum = assert_eq(ndum,size(q),size(LHM,1)-2,'setLHM 2')
ndum = assert_eq(ndum,size(R),size(diagk),'setLHM 3')
ndum = assert_eq(ndum,size(LHM,2)-2,'setLHM 4')

LHM=0.0_dp

J = DiagMul(R,J,ndum) !! Destroys J, but we only need RJ

LHM(1:size(v),1:size(v)) = UnitMatrix(size(v)) - J 

LHM(size(LHM,1)-1,1:size(q)) = MatMul(q,J) 
LHM(size(LHM,1),1:size(q)) = q

LHM(1:size(v),size(LHM,1)-1) = -R*v
LHM(1:size(v),size(LHM,1)) = -R*(diagk*f)

LHM(size(LHM,1)-1,size(LHM,1)-1) = -dot_product(q,LHM(1:size(v),size(LHM,1)-1))
LHM(size(LHM,1)-1,size(LHM,1)) = -dot_product(q,LHM(1:size(v),size(LHM,1)))


end subroutine