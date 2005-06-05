Subroutine eulerJ(x,y,J,Jout,h,MatVar)

USE nrtype

IMPLICIT NONE

REAL(DP), DIMENSION(:), INTENT(IN) :: y
REAL(DP), INTENT(IN) :: x,h
REAL(DP), DIMENSION(:,:), INTENT(IN) :: J
REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Jout
interface
	function MatVar(x,y)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(size(y),size(y)) :: MatVar
	end function MatVar
end interface

Jout=J+h*matmul(MatVar(x,y),J)

end subroutine