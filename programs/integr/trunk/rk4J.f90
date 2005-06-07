Subroutine rk2J(x,y,dydx,h,yout,J,dJds,Jout,MatVar,derivs)

USE nrtype ; USE ifc_integr, ONLY: EulerJ
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
	function MatVar(x,y)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(size(y),size(y)) :: MatVar
	end function MatVar
end interface


REAL(DP), DIMENSION(size(y)) :: k1,k2,v
real(dp), dimension(size(J,1),size(J,2)) ::kJ1, kJ2,kJ3,kJ4,Jdum,dJdsdum 
integer(i4b) :: ndum

ndum=assert_eq(size(y),size(dydx),size(yout),'rk2J')


k1 = h*dydx ! First step
kJ1 = h*dJds
call derivs(x + h/2.0_dp , y + k1/2.0_dp, v)  ! second step
k2 = h*v
call eulerJ(x,y(1:size(y)-1),J,Jdum,h/2.0_dp, MatVar)
call derivsJ(x + h/2.0_dp , y + k1/2.0_dp, Jdum, dJdsdum, MatVar)
kJ2=h*dJdsdum
call derivs(x + h/2.0_dp , y + k2/2.0_dp, v, kappa) ! third step
k3 = h*v
call eulerJ(x + h/2.0_dp , y + k2/2.0_dp,J,Jdum,h/2.0_dp, MatVar)
call derivsJ(x + h/2.0_dp , y + k2/2.0_dp, Jdum, dJdsdum, MatVar)
kJ3=h*dJdsdum
call derivs(x + h, y + k3,v)
k4=h*v 
call eulerJ(x + h, y + k3,J,Jdum,h, MatVar)
call derivsJ(x + h, y + k3, Jdum, dJdsdum, MatVar)
kJ4=h*dJdsdum

yout = y + k1/6.0_dp+ k2/3.0_dp+ k3/3.0_dp+ k4/6.0_dp	! Accumulate increments
Jout = J + kJ1/6.0_dp+ kJ2/3.0_dp+ kJ3/3.0_dp+ kJ4/6.0_dp	


end subroutine