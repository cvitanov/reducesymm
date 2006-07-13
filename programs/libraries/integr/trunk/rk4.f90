Subroutine rk4(x,y,dydx,h,yout,derivs)

! From Numerical Recipes Book.


USE nrtype; USE nrutil, ONLY : assert_eq

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

!Given values for the N variables y and their derivatives dydx known at x, use the fourthorder
!Runge-Kutta method to advance the solution over an interval h and return the incremented
!variables as yout, which need not be a distinct array from y. y, dydx and yout
!are all of length N. The user supplies the subroutine derivs(x,y,dydx), which returns
!derivatives dydx at x.

INTEGER(I4B) :: ndum
REAL(dp) :: h6,hh,xh
REAL(dp), DIMENSION(size(y)) :: dym,dyt,yt

ndum=assert_eq(size(y),size(dydx),size(yout),'rk4')

hh=h*0.5_dp
h6=h/6.0_dp
xh=x+hh
yt=y+hh*dydx !First step.

call derivs(xh,yt,dyt) !Second step.
yt=y+hh*dyt
call derivs(xh,yt,dym) !Third step.
yt=y+h*dym
dym=dyt+dym
call derivs(x+h,yt,dyt) !Fourth step.
yout=y+h6*(dydx+dyt+2.0_dp*dym) !Accumulate increments with proper weights.



END SUBROUTINE rk4
