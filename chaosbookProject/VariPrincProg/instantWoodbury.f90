SUBROUTINE instantWoodbury(Usaved,V,Hinv,b,x)

USE interfcs
USE nrtype

IMPLICIT NONE

REAL(DP), INTENT(INOUT) ::  Usaved(:,:)
REAL(DP), INTENT(IN) :: V(:,:)
REAL(DP), INTENT(INOUT) :: Hinv(:,:)
REAL(DP), INTENT(IN) :: b(:)
REAL(DP), INTENT(OUT) :: x(:)



! Apply Woodbury formula without storage of the inverse



!Calculate the solution to the initial problem. b is y in nr book
x=b-MatMul(Usaved,MatMul(Hinv,MatMul(Transpose(V),b) ) )



END SUBROUTINE