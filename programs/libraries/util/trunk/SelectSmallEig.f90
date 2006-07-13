logical function SelectSmallEig(wR_j,wI_j)

use nrtype
implicit none

real(dp), INTENT(IN) :: wR_j,wI_j
!
!
real(dp) :: M

M=100.0
if ( abs(wR_j+ii*wI_j) < M ) then
	SelectSmallEig=.true.
else
	SelectSmallEig=.false.
endif

end function