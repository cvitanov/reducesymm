program roesslerP

use nrtype
use ifc_integr, only: integrP

implicit none

integer(i4b) :: d=3, nsteps=1000, nstepsP=2000, nInters=100, sect=1
real(dp), allocatable :: yP(:,:), yi(:)
real(dp) ::  Delta_x=0.1_dp, qfP=0.0_dp
INTERFACE
	SUBROUTINE roesslerField(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp), DIMENSION(:), INTENT(IN) :: y
		REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE roesslerField
END INTERFACE


integer(i4b) :: i

allocate(yP(nInters,d+1), yi(d+1))

yi(1)=-0.1_dp
yi(2)= 0.0_dp
yi(3)= 0.0_dp
yi(d+1)= 0.0_dp

call integrP(yi,Delta_x,qfP,yP,nsteps,nstepsP,nInters,sect,roesslerField)

open(9,file="roesslerP.dat")
do i=1,nInters
	write (9,"4F12.7") yP(i,:)
end do
close(9)


end program