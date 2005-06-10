program newtonSearch

use nrtype; use parameters
use ifc_newton
use ifc_integr
use ifc_util
use F95_LAPACK, only: LA_GEEV

implicit none

real(dp) :: T, yi(d+1), a(d), xi, J(d,d)
integer(i4b) :: conv=0
real(dp), dimension(d) :: WI, WR


INTERFACE
	SUBROUTINE roesslerField(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE roesslerField
END INTERFACE
interface
	function roesslerVar(x,y)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(size(y),size(y)) :: roesslerVar
	end function 
end interface


T=Tguess

call init_a(a)

open (9, file='guess.dat')
read(9,*) yi(1:d)
close(9)

yi(d+1)=0.0_dp

call NewtonPO(yi,T,a,tol,maxIter,nstepsN,roesslerField,roesslerVar,conv)

print *,conv

xi=0.0_dp

!call integrP(yi,Delta_x,qfP,yP,nsteps,nstepsP,nInters,sect,derivs)

J=UnitMatrix(d)

call rk4Jdriver(xi,yi,T,nsteps,yi,J,J,roesslerVar,roesslerField)

call LA_GEEV( J, WR, WI)

Print *, WR
Print *,WI

Print '(E1.3)',1.0_dp-WR(d)

open (9, file='cycle.dat')
write(9,format_label) yi(1:d)
close(9)

open (9, file='period.dat')
write(9,format_label) T
close(9)

end program
