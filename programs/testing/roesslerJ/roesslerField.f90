SUBROUTINE roesslerField(x,y,dydx)

USE nrtype

IMPLICIT NONE

REAL(DP), INTENT(IN) :: x
REAL(DP), DIMENSION(:), INTENT(IN) :: y
REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx

REAL(DP) :: a=0.2,b=0.2,c=5.7


dydx(1)= ( -y(2)-y(3) )
dydx(2)= ( y(1) + a*y(2) )
dydx(3)= ( b + y(3)*( y(1) - c ) )
dydx(4)= 1.0_dp

END SUBROUTINE