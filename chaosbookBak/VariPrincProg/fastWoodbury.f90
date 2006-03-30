SUBROUTINE fastWoodbury(LHM,Usaved,V,Hinv,N,d,m1,m2,b,al,indx,x)

USE interfcs
USE nrtype

IMPLICIT NONE

REAL(DP), INTENT(IN) :: LHM(:,:)
REAL(DP), INTENT(INOUT) ::  Usaved(:,:)
REAL(DP), INTENT(IN) :: V(:,:)
REAL(DP), INTENT(IN) :: Hinv(:,:)
INTEGER(I4B), INTENT(IN) :: N,m1,m2,d
REAL(DP), INTENT(IN) :: al(:,:)
INTEGER(I4B), INTENT(IN) :: indx(:)
REAL(DP), INTENT(IN) :: b(:)
REAL(DP), INTENT(OUT) :: x(:)



! Apply Woodbury formula without storage of the inverse. "fast" version


INTEGER(I4B) ::  i


! The rest of the rhs columns being the same, only solve the auxiliary problem
call banbks(LHM,m1,m2,al,indx,Usaved(:,4*d+1)) ! U is both input and output (i.e. Z in nr book)

!Calculate the solution to the initial problem. b is y in nr book
x=b-MatMul(Usaved,MatMul(Hinv,MatMul(Transpose(V),b) ) )


END SUBROUTINE