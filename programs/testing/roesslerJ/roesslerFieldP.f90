SUBROUTINE roesslerFieldP(x,y,dydx,kappa)

USE nrtype

IMPLICIT NONE

REAL(DP), INTENT(IN) :: x
REAL(DP), DIMENSION(:), INTENT(IN) :: y
REAL(DP), INTENT(IN) :: kappa
REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx

REAL(DP) :: a=0.2,b=0.2,c=5.7


dydx(1)= kappa*( -y(2)-y(3) )
dydx(2)= kappa*( y(1) + a*y(2) )
dydx(3)= kappa*( b + y(3)*( y(1) - c ) )
dydx(4)= kappa

END SUBROUTINE