logical function SelectSmallEig_r(wR_j,wI_j)

use nrtype
implicit none

real(dp), INTENT(IN) :: wR_j,wI_j
!
!
real(dp) :: M

M=100.0
if ( abs(wR_j+ii*wI_j) < M ) then
	SelectSmallEig_r=.true.
else
	SelectSmallEig_r=.false.
endif

end function SelectSmallEig_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

logical function SelectSmallEig_c(w_j)

use nrtype
implicit none

complex(dpc), INTENT(IN) :: w_j
!
!
real(dp) :: M

M=10.0
if ( abs(w_j) < M ) then
	SelectSmallEig_c=.true.
else
	SelectSmallEig_c=.false.
endif

end function SelectSmallEig_c
