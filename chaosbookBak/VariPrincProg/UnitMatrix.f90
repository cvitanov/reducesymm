FUNCTION UnitMatrix(N)

USE nrtype

IMPLICIT NONE

INTEGER(I4B) :: N
REAL(DP) :: UnitMatrix(N,N)

!Returns the NxN identity matrix

INTEGER(I4B) :: i

UnitMatrix=0

do i=1,size(UnitMatrix,1)
  UnitMatrix(i,i)=1
end do

END FUNCTION