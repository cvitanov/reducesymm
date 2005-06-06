SUBROUTINE roesslerField(x,y,dydx)

USE nrtype
USE parameters

IMPLICIT NONE

REAL(DP), INTENT(IN) :: x
REAL(DP), DIMENSION(:), INTENT(IN) :: y
REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx

dydx(1)= ( -y(2)-y(3) )
dydx(2)= ( y(1) + alpha*y(2) )
dydx(3)= ( beta + y(3)*( y(1) - gamma ) )
dydx(4)= 1.0_dp

END SUBROUTINE