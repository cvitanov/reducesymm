SUBROUTINE roesslerFieldP(x,y,dydx,kappa)

USE nrtype
USE parameters

IMPLICIT NONE

REAL(DP), INTENT(IN) :: x
REAL(DP), DIMENSION(:), INTENT(IN) :: y
REAL(DP), INTENT(IN) :: kappa
REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx

dydx(1)= kappa*( -y(2)-y(3) )
dydx(2)= kappa*( y(1) + alpha*y(2) )
dydx(3)= kappa*( beta + y(3)*( y(1) - gamma ) )
dydx(4)= kappa

END SUBROUTINE