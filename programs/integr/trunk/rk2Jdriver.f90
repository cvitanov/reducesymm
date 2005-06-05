SUBROUTINE rk2Jdriver(xi,yi,xf,nsteps,y,Ji,Jout,MatVar,derivs)

USE nrtype ; USE ifc_integr, ONLY: rk2J 
USE nrutil, ONLY : assert_eq

implicit none

REAL(DP), INTENT(IN) :: xi,xf
REAL(DP), DIMENSION(:), INTENT(IN) :: yi
REAL(DP), DIMENSION(:,:), INTENT(OUT) :: y
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

REAL(DP) :: h, v(size(yi))
INTEGER(I4B) :: i,ndum, mdum

ndum=assert_eq(size(y,1),nsteps+1,'rk2Jdriver: y dim1')
mdum=assert_eq(size(y,2),size(yi),'rk2Jdriver: y dim2')

h=(xf-xi)/Real(nsteps,dp)

y(1,:)=yi

DO i=2,nsteps+1
	call derivs(y(i-1,size(y,2)),y(i-1,:),v)
	call rk2J(y(i-1,size(y,2)),y(i-1,:),v,h,y(i,:),Ji,Jout,MatVar,derivs)
END DO





END SUBROUTINE