
! Fourth order Runge-Kutta integrator. Needs the module nrtype from Numerical Recipes,
! that defines data types. Time is treated as an extra variable of the system 
! so that the program needs fewer modifications to construct Poincare sections.
! Here a version to integrate a system of ODE's defined in a complex space is given.
! Uses KS equations as example. 


PROGRAM integr

use parameters
USE nrtype
USE ifc_integr

IMPLICIT NONE

REAL(dp), DIMENSION(:,:), ALLOCATABLE :: y
REAL(dp) :: xi=0.0_dp, xf=10.0_dp
INTEGER(I4B) :: nsteps=100000, nrep =10,  d=16,i


INTERFACE
	SUBROUTINE KSfield(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp), DIMENSION(:), INTENT(IN) :: y
		REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE KSfield
END INTERFACE

allocate(y(nsteps+1,2*d+1))

y=0.0_dp

OPEN(10, file='ic.dat')

	read(10,*) y(1,d+1:2*d) 

CLOSE(10)

y(1,2*d+1)=xi

y(1,1)=1.0_dp
y(1,2)=1.0_dp
y(1,4)=1.0_dp
y(1,6)=1.0_dp
y(1,9)=1.0_dp

do i=1,nrep 
	call rk4driver(xi,y(1,:),xf,nsteps,y,KSfield)
	y(1,:)=y(size(y,1),:) 
end do

OPEN(10, file='ksR.dat')

DO i=1, size(y,1)
	WRITE(10,format_label) y(i,1:d) 
END DO

CLOSE(10) 

OPEN(10, file='ksI.dat')

DO i=1, size(y,1)
	WRITE(10,format_label) y(i,d+1:2*d) 
END DO

CLOSE(10) 

OPEN(10, file='kst.dat')

DO i=1, size(y,1)
	WRITE(10,format_label) y(i,2*d+1) 
END DO

CLOSE(10) 


END PROGRAM
