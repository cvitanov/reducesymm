subroutine rk4J_nr(x,y,dydx,h,yout,J,dJds,Jout,derivs,derivsJ)

USE nrtype !; USE ifc_integr, ONLY: derivsJ
USE nrutil, ONLY: assert_eq

IMPLICIT NONE

REAL(DP), DIMENSION(:), INTENT(IN) :: y,dydx
REAL(DP), INTENT(IN) :: x,h
REAL(DP), DIMENSION(:), INTENT(OUT) :: yout
REAL(DP), DIMENSION(:,:), INTENT(IN) :: J, dJds
REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Jout
INTERFACE
	SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE derivs
END INTERFACE
interface
	SUBROUTINE derivsJ(s,J,dJds,a)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: s
		REAL(dp), DIMENSION(:,:), INTENT(IN) :: J
		REAL(dp), DIMENSION(:,:), INTENT(OUT) :: dJds
		REAL(dp), DIMENSION(:), INTENT(IN) :: a
	end subroutine
end interface
INTEGER(I4B) :: ndum,k
REAL(DP) :: h6,hh,xh
REAL(DP), DIMENSION(size(y)) :: dym,dyt,yt
REAL(DP), DIMENSION(size(J,1),SIZE(J,1)) :: dJm,dJt,Jt

ndum=assert_eq(size(y),size(dydx),size(yout),'rk4J')

hh=h*0.5_dp
h6=h/6.0_dp

xh=x+hh
yt=y+hh*dydx !First step.
Jt=J+hh*dJdx
call derivs(xh,yt,dyt)  ! second step
call derivsJ(xh,Jt,dJt,yt)
yt=y+hh*dyt
Jt=J+hh*dJt
call derivs(xh,yt,dym) !Third step.
call derivsJ(xh,Jt,dJm,yt)
yt=y+h*dym
Jt=J+h*dJm
dym=dyt+dym
dJm=dJt+dJm
call derivs(x+h,yt,dyt) !Fourth step.
call derivsJ(x+h,Jt,dJt,yt)

yout=y+h6*(dydx+dyt+2.0_dp*dym) !Accumulate increments with proper weights.
Jout=J+h6*(dJdx+dJt+2.0_dp*dJm)

end subroutine