! Exercise 2.6

! Fourth order Runge-Kutta integrator. Needs the module nrtype from Numerical Recipes,
! that defines data types. Time is treated as an extra variable of the system 
! so that the program needs fewer modifications to construct Poincare sections.



PROGRAM integr


USE nrtype; USE ifc_integr

IMPLICIT NONE

REAL(dp), DIMENSION(:,:), ALLOCATABLE :: y
REAL(dp) :: xi=0, xf=10
INTEGER(I4B) :: nsteps=200, dmn=3,i

INTERFACE
	SUBROUTINE oscillator(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp), DIMENSION(:), INTENT(IN) :: y
		REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE oscillator
END INTERFACE



allocate(y(nsteps+1,dmn))

y=0

y(1,1)=xi
y(1,2)=0
y(1,3)=1

call rk4driver(xi,y(1,:),xf,nsteps,y,oscillator) 

OPEN(10, file='oscillator.dat')

DO i=1, size(y,1)
	WRITE(10,*) y(i,:) 
END DO

CLOSE(10) 


END PROGRAM