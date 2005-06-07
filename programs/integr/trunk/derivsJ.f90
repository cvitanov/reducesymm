Subroutine derivsJ(x,y,J,dJds,MatVar)

USE nrtype

IMPLICIT NONE

REAL(DP), DIMENSION(:), INTENT(IN) :: y
REAL(DP), INTENT(IN) :: x
REAL(DP), DIMENSION(:,:), INTENT(IN) :: J
REAL(DP), DIMENSION(:,:), INTENT(OUT) ::dJds
interface
	function MatVar(x,y)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(size(y),size(y)) :: MatVar
	end function MatVar
end interface

dJds = matmul(MatVar(x,y),J)

end subroutine