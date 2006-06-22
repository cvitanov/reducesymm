Subroutine derivsJ(x,y,J,dJdx,MatVar)

USE nrtype
USE nrutil, ONLY : assert_eq

IMPLICIT NONE

REAL(DP), DIMENSION(:), INTENT(IN) :: y
REAL(DP), INTENT(IN) :: x
REAL(DP), DIMENSION(:,:), INTENT(IN) :: J
REAL(DP), DIMENSION(:,:), INTENT(OUT) ::dJdx
interface
	function MatVar(x,y)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(size(y),size(y)) :: MatVar
	end function MatVar
end interface

integer(i4b) ndum

ndum=assert_eq(size(y),size(J,1),size(dJdx,1),'derivsJ')

dJdx = matmul(MatVar(x,y),J)

end subroutine