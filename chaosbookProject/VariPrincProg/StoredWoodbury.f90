! author: Evangelos Siminos, School of Physics, Georgia Institute of Technology

SUBROUTINE StoredWoodbury(aInv,U,V)


USE nrtype; USE interfcs, ONLY: LUinv


IMPLICIT NONE

REAL(DP), INTENT(INOUT) :: aInv(:,:)
REAL(DP), INTENT(IN) :: U(:,:), V(:,:)

! This routine uses Woodbury formula to compute the correction to the inverse aInv 
! of an N x N matrix a due to a perturbation that has the form U.Transpose(V),  
! where U and V are N x P matrices with P<<N. The updated inverse is returned in
! aInv.


INTEGER(I4B) :: N, P,i,j
REAL(DP), ALLOCATABLE :: Cor1(:,:), UnitMatrix(:,:)	

N=size(aInv,1)
P=size(U,2)											

allocate(Cor1(P,P))
allocate(UnitMatrix(P,P))

UnitMatrix=0
forall (i=1:P, j=1:P, i==j) UnitMatrix(i,j)=1


!Calculate central term
Cor1 = UnitMatrix + MatMul( transpose(V), MatMul(aInv, U) )


!and its inverse needed by W.F.:

!Print *,"r",size(Cor1,1)

if (size(Cor1,1)==1) then ! if Cor1 is scalar then we have the Sherman-Morrison case
	Cor1(1,1)=1/Cor1(1,1)	
else					! a true matrix neads to be inverted
	call LUinv(Cor1)
end if

!Apply Woodbury Formula
aInv = aInv - MatMul(aInv, MatMul(U,MatMul(Cor1, MatMul(Transpose(V),aInv) ) ) )

!Print *, aInv


END SUBROUTINE