MODULE ifc_integr

! Contains all the subroutine declarations.

INTERFACE
	SUBROUTINE rk4(y,dydx,x,h,yout,derivs)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), DIMENSION(:), INTENT(IN) :: y,dydx
		REAL(dp), INTENT(IN) :: x,h
		REAL(dp), DIMENSION(:), INTENT(OUT) :: yout
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
				USE nrtype
				IMPLICIT NONE
				REAL(dp), INTENT(IN) :: x
				REAL(dp), DIMENSION(:), INTENT(IN) :: y
				REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx	
			END SUBROUTINE derivs
		END INTERFACE
	END SUBROUTINE
END INTERFACE

INTERFACE
	SUBROUTINE rk4driver(xi,yi,xf,nsteps,y,derivs)
		USE nrtype
		IMPLICIT none
		REAL(dp), INTENT(IN) :: xi,xf
		REAL(dp), DIMENSION(:), INTENT(IN) :: yi
		REAL(dp), DIMENSION(:,:), INTENT(OUT) :: y
		INTEGER(I4B), INTENT(IN) :: nsteps
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
				USE nrtype
				IMPLICIT NONE
				REAL(dp), INTENT(IN) :: x
				REAL(dp), DIMENSION(:), INTENT(IN) :: y
				REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx	
			END SUBROUTINE derivs
		END INTERFACE
	END SUBROUTINE
END INTERFACE


END MODULE