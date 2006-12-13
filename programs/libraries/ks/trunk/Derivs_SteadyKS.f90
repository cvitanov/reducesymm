SUBROUTINE Derivs_SteadyKS(x,y,dydx)

USE nrtype
use ks
use nrutil

IMPLICIT NONE
REAL(dp), INTENT(IN) :: x
REAL(dp), DIMENSION(:), INTENT(IN) :: y
REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx
!
integer(i4b) :: dim

dim=assert_eq(size(y),size(dydx),3,'Derivs_SteadyKS')

dydx(1)=y(2)
dydx(2)=y(3)
dydx(3)=y(1)**2-y(2)-c

END SUBROUTINE 
