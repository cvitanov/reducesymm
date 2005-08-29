! Exercise 2.6

! Fourth order Runge-Kutta integrator. Needs the module nrtype from Numerical Recipes,
! that defines data types. Time is treated as an extra variable of the system 
! so that the program needs fewer modifications to construct Poincare sections.



PROGRAM integr


USE nrtype; USE ifc_integr; USE ode_path

IMPLICIT NONE

REAL(dp), DIMENSION(:), ALLOCATABLE :: y
REAL(dp) :: xi=0, xf=10
INTEGER(I4B) :: dmn=3,i
real(dp) :: eps=1e-10_dp, h1=1e-3_dp, hmin=1e-14_dp

INTERFACE
	SUBROUTINE oscillator(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp), DIMENSION(:), INTENT(IN) :: y
		REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE oscillator
END INTERFACE

save_steps=.true.
dxsav = 1e-6_dp


allocate(y(dmn))

y=0

y(1)=xi
y(2)=0
y(3)=1

call odeint(y,xi,xf,eps,h1,hmin,oscillator,rkqs) 

OPEN(10, file='oscillator.dat')

DO i=1, size(yp,2)
	WRITE(10,*) yp(:,i) 
END DO

CLOSE(10) 


END PROGRAM