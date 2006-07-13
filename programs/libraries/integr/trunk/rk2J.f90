Subroutine rk2J(x,y,dydx,h,yout,J,Jout,MatVar,derivs)

USE nrtype ; USE ifc_integr, ONLY: EulerJ
USE nrutil, ONLY: assert_eq

IMPLICIT NONE

REAL(DP), DIMENSION(:), INTENT(IN) :: y,dydx
REAL(DP), INTENT(IN) :: x,h
REAL(DP), DIMENSION(:), INTENT(OUT) :: yout
REAL(DP), DIMENSION(:,:), INTENT(IN) :: J
REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Jout
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


REAL(DP), DIMENSION(size(y)) :: k1,k2,v
real(dp), dimension(size(J,1),size(J,2)) ::kJ1, kJ2
integer(i4b) :: ndum

ndum=assert_eq(size(y),size(dydx),size(yout),'rk2J')

kJ1=0.0_dp
kJ2=0.0_dp

k1 = h*dydx
call eulerJ(x,y(1:size(y)-1),J,kJ1,h/2.0_dp, MatVar)
call derivs(x + h/2.0_dp , y + k1/2.0_dp, v)
kJ2 = h*matmul(MatVar(x+h/2.0_dp,y(1:size(y)-1) + k1(1:size(y)-1)/2.0_dp),kJ1)
k2 = h*v
yout = y+k2
Jout = J+kJ2


end subroutine