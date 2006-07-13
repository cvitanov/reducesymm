logical function SelectLargeEig(wR_j,wI_j)

use nrtype
implicit none

real(dp), INTENT(IN) :: wR_j,wI_j
!
!
REAL(dp) :: eps

eps=epsilon(1.0_dp)!real(ibeta,dp)**machep
if ( abs(wR_j+ii*wI_j)>eps ) then
	SelectLargeEig=.true.
else
	SelectLargeEig=.false.
endif

end function