SUBROUTINE rk4driver(xi,yi,xf,nsteps,y,derivs)

! Driver that simply repeats the integration steps of rk4 in order
! to cover a range xi to xf of the independed variable.


USE nrtype ; USE ifc_integr

IMPLICIT none

REAL(dp), INTENT(IN) :: xi,xf
REAL(dp), DIMENSION(:), INTENT(IN) :: yi
REAL(dp), DIMENSION(:,:), INTENT(OUT) :: y
INTEGER(I4B), INTENT(IN) :: nsteps
INTERFACE
	SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp), DIMENSION(:), INTENT(IN) :: y
		REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE derivs
END INTERFACE
!
!

REAL(dp) ::h, v(size(yi))
INTEGER(I4B) ::i

h = (xf-xi)/nsteps

y(1,:)=yi

DO i=2,nsteps
	call derivs(y(i-1,1),y(i-1,:),v)
	CALL rk4(y(i-1,:),v,y(i-1,1),h,y(i,:),derivs)
END DO





END SUBROUTINE