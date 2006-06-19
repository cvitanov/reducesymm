SUBROUTINE Woodbury(LHM,U,V,N,d,m1,m2,b,al,indx,Ho,Hinv,deltaA)

USE interfcs
USE nrtype

IMPLICIT NONE

REAL(DP), INTENT(IN) :: LHM(:,:)
REAL(DP), INTENT(INOUT) ::  U(:,:)
REAL(DP), INTENT(IN) :: V(:,:)
INTEGER(I4B), INTENT(IN) :: N,m1,m2,d
REAL(DP), INTENT(IN) :: al(:,:)
INTEGER(I4B), INTENT(IN) :: indx(:)
REAL(DP), INTENT(IN) :: b(:)
REAL(DP), INTENT(OUT) :: Ho(:,:)
REAL(DP), INTENT(OUT) :: Hinv(:,:)
REAL(DP), INTENT(OUT) :: deltaA(:)




! Apply Woodbury formula without storage of the inverse


INTEGER(I4B) ::  i


! Solve the 4*d+1 auxiliary problems
do i=1,4*d+1
	call banbks(LHM,m1,m2,al,indx,U(:,i)) ! U is both input and output (i.e. Z in nr book)
end do
 
! Find the auxiliary matrix H (notice that now U is Z in nr book)
Ho=UnitMatrix(4*d+1)+MatMul(Transpose(V),U)
Print *, "before wood inversion"
Hinv=Ho
call LUinv(Hinv)
Print *, "after wood inversion"

!Calculate the solution to the initial problem. b is y in nr book
deltaA=b-MatMul(U,MatMul(Hinv,MatMul(Transpose(V),b) ) )


END SUBROUTINE