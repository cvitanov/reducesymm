subroutine setLHM(f,J,v,q,LHM)

use nrtype 
use nrutil, only: assert_eq
use ifc_util, only: UnitMatrix

implicit none

complex(dpc), dimension(:,:), intent(in) :: J
complex(dpc), dimension(:), intent(in) ::  f, v, q
complex(dpc), dimension(:,:), intent(out) :: LHM 
!!
integer(i4b) :: ndum

ndum = assert_eq(size(v),size(J,1),size(f),'setLHM 1')
ndum = assert_eq(ndum,size(q)/2-1,size(LHM,1)/2-1,'setLHM 2')
ndum = assert_eq(ndum,size(LHM,2)/2-1,'setLHM 3')


LHM=0.0_dp

LHM(1:size(v),1:size(v)) = UnitMatrix(size(v)) - J 

LHM(size(J,1)+1,1:size(q)) = q

LHM(1:size(v),size(J,2)+1) = 0.5_dp*v

LHM(1:size(v),size(J,2)+2) = 0.5_dp*v

!!! The imaginary part of the period should remain zero
LHM(size(J,1)+2,size(J,2)+1) = (1.0,0.0)
LHM(size(J,1)+2,size(J,2)+2) = -(1.0,0.0)
!!!

LHM(size(v)+3:size(LHM,1),size(J,2)+1) = 0.5_dp*conjg(v)

LHM(size(v)+3:size(LHM,1),size(J,2)+2) = 0.5_dp*conjg(v)

LHM(size(v)+3:size(LHM,1),size(v)+3:size(LHM,2)) = UnitMatrix(size(v)) - conjg(J) 



end subroutine