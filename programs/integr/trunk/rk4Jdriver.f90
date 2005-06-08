SUBROUTINE rk4Jdriver(xi,yi,xf,nsteps,y,Ji,Jout,MatVar,derivs)

USE nrtype ; USE ifc_integr, ONLY: rk4J, derivsJ 
USE nrutil, ONLY : assert_eq

implicit none

REAL(DP), INTENT(IN) :: xi,xf
REAL(DP), DIMENSION(:), INTENT(IN) :: yi
REAL(DP), DIMENSION(:), INTENT(OUT) :: y
INTEGER(I4B), INTENT(IN) :: nsteps
REAL(DP), DIMENSION(:,:), INTENT(IN) :: Ji
real(dp), dimension(:,:), intent(out) :: Jout
INTERFACE
	SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE derivs	
END INTERFACE
interface
	function MatVar(x,y)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(size(y),size(y)) :: MatVar
	end function MatVar
end interface
!
!

REAL(DP) :: h, v(size(yi)), dJds(size(Ji,1),size(Ji,2))
INTEGER(I4B) :: i, mdum

mdum=assert_eq(size(y),size(yi),'rk2Jdriver: y dim2')

h=(xf-xi)/Real(nsteps,dp)
y(:)=yi

DO i=2,nsteps+1
	call derivs(y(size(y)),y,v)
	call derivsJ(y(size(y)),y(1:size(y)-1), Ji, dJds, MatVar)
	call rk4J(y(size(y)),y,v,h,y,Ji,dJds,Jout,MatVar,derivs)
END DO



END SUBROUTINE