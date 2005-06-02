SUBROUTINE rk4Pdriver(xi,yi,xf,nsteps,y,derivs,p,sect)

USE nrtype ; USE ifc_integr, ONLY : rk4P ; USE nrutil, ONLY : assert_eq

IMPLICIT none

REAL(DP), INTENT(IN) :: xi,xf
REAL(DP), DIMENSION(:), INTENT(IN) :: yi
REAL(DP), DIMENSION(:,:), INTENT(OUT) :: y
INTEGER(I4B), INTENT(IN) :: nsteps, p, sect
INTERFACE
	SUBROUTINE derivs(x,y,dydx,kappa)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x,kappa
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE derivs
END INTERFACE
!
!

REAL(DP) ::h, v(size(yi)), kappa
INTEGER(I4B) ::i,idum, ndum, mdum

ndum=assert_eq(size(y,1),nsteps+1,'rk4Pdriver: y dim1')
mdum=assert_eq(size(y,2),size(yi),'rk4Pdriver: y dim2')

h=(xf-xi)/Real(nsteps,dp)

y(1,:) = yi

kappa = 1.0_dp


if ( p == 1 ) then
	DO i=2,nsteps+1
		kappa=1.0_dp	
		call derivs(y(i-1,size(y,2)),y(i-1,:),v, kappa)
		kappa=1.0_dp/v(sect)
		call derivs(y(i-1,size(y,2)),y(i-1,:),v, kappa)
		CALL rk4P(y(i-1,:),v,y(i-1,size(y,2)),h,y(i,:),derivs,p,sect)
	END DO
else
	DO i=2,nsteps+1
		call derivs(y(i-1,size(y,2)),y(i-1,:),v, kappa)
		CALL rk4P(y(i-1,:),v,y(i-1,size(y,2)),h,y(i,:),derivs,p,sect)
	END DO
end if




END SUBROUTINE