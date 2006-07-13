FUNCTION UnitMatrix(N)

USE nrtype

IMPLICIT NONE

INTEGER(I4B) :: N
REAL(DP) :: UnitMatrix(N,N)

!Returns the NxN identity matrix

INTEGER(I4B) :: i

UnitMatrix=0.0_dp

do i=1,size(UnitMatrix,1)
  UnitMatrix(i,i)=1.0_dp
end do

END FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
