MODULE intfaces

INTERFACE
	SUBROUTINE roesslerField(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE roesslerField
END INTERFACE

INTERFACE
	SUBROUTINE roesslerFieldP(x,y,dydx,kappa)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x,kappa
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE roesslerFieldP
END INTERFACE

interface
	function roesslerVar(x,y)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(size(y),size(y)) :: MatVar
	end function roesslerVar
end interface

INTERFACE
	FUNCTION UnitMatrix(N)
		USE nrtype
		IMPLICIT NONE
		INTEGER(I4B) :: N
		REAL(DP) :: UnitMatrix(N,N)
	END FUNCTION
END INTERFACE

END MODULE