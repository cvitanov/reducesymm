SUBROUTINE oscillator(x,y,dydx)

USE nrtype

IMPLICIT NONE

REAL(dp), INTENT(IN) :: x
REAL(dp), DIMENSION(:), INTENT(IN):: y
REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx

dydx(1)=1
dydx(2)=y(3)
dydx(3)=-y(2)


END SUBROUTINE