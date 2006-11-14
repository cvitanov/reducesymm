logical function SelectLargeEig_r(wR_j,wI_j)

use nrtype
implicit none

real(dp), INTENT(IN) :: wR_j,wI_j
!
!
REAL(dp) :: eps

eps=epsilon(1.0_dp)!real(ibeta,dp)**machep
if ( abs(wR_j+ii*wI_j)>eps ) then
	SelectLargeEig_r=.true.
else
	SelectLargeEig_r=.false.
endif

end function SelectLargeEig_r

logical function SelectLargeEig_c(w_j)

use nrtype
implicit none

complex(dpc), INTENT(IN) :: w_j
!
!
REAL(dp) :: eps

eps=epsilon(1.0_dp)!real(ibeta,dp)**machep
if ( abs(w_j)>eps ) then
	SelectLargeEig_c=.true.
else
	SelectLargeEig_c=.false.
endif

end function SelectLargeEig_c