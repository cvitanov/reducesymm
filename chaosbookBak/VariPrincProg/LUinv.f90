SUBROUTINE LUinv(a)


USE nrtype
USE interfcs, ONLY: ludcmp,lubksb

IMPLICIT NONE

REAL(DP), INTENT(INOUT) :: a(:,:)


! Uses LU-decomposition routines ludcmp and lubksb from Numerical Recipes to invert a 
! matrix a. The inverse is returned in place of a which is destroyed (can be modified). 


REAL(DP) :: d
INTEGER(I4B)  :: N,i,j 
INTEGER(I4B), allocatable  :: indx(:)
REAL(DP), allocatable, DIMENSION(:,:) :: aInv

N=size(a,1)

allocate( indx(N), aInv(N,N) )


!LU-decompose it
call ludcmp(a,indx,d)
!and invert it column by column
	! Form the unit matrix
	aInv=0
	forall (i=1:N, j=1:N, i==j) aInv(i,j)=1
do j=1,N
	call lubksb(a,indx,aInv(:,j))
end do

a=aInv

END SUBROUTINE