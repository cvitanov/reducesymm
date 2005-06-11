program newtonSearch

use nrtype; use parameters
use ifc_newton, only: NewtonPO
use ifc_integr, only: rk4Jdriver, integrP
use ifc_util, only: UnitMatrix
use F95_LAPACK, only: LA_GEEV

implicit none

real(dp) :: T, yi(d+1), y(d+1), a(d), xi, J(d,d), qfP, yP(nInters,d+1)
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
INTERFACE
	SUBROUTINE roesslerFieldP(x,y,dydx,kappa)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x, kappa
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE roesslerFieldP
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

call NewtonPO(yi,T,a,tol,maxIter,nstepsN,roesslerVar,roesslerField,conv,J)

print *,"converged?", conv

call LA_GEEV( J, WR, WI)

print *,"Eigenvalues of J, as returned from NewtonPO"
Print *,"R",WR
Print *,"I",WI



xi=0.0_dp

!call integrP(yi,Delta_x,qfP,yP,nsteps,nstepsP,nInters,sect,derivs)

J=UnitMatrix(d)

call rk4Jdriver(xi,yi,T,nsteps,y,J,J,roesslerVar,roesslerField)

Print *, "distance from initial point after integration:"
print *, y(1:d)-yi(1:d)

call LA_GEEV( J, WR, WI)

print *,"Eigenvalues of J"
Print *,"R",WR
Print *,"I",WI

Print *,"marginal eigenvalue-1",1.0_dp-WR(d)

qfP = yi(sect)

call integrP(yi,Delta_x,qfP,yP,nstepsN,nstepsP,nInters,sect,direction,roesslerFieldP)

print *, "Distance from initial point after integrating with integrP"
print *, yP(1,1:d)-yi(1:d)

J=UnitMatrix(d)

T = yP(d+1)

call rk4Jdriver(xi,yi,T,nsteps,y,J,J,roesslerVar,roesslerField)

print *,"Eigenvalues of J"
Print *,"R",WR
Print *,"I",WI

Print *,"marginal eigenvalue-1",1.0_dp-WR(d)


open (9, file='cycle.dat')
write(9,format_label3) yi(1:d)
close(9)

open (9, file='period.dat')
write(9,'(F20.17)') T
close(9)

end program
