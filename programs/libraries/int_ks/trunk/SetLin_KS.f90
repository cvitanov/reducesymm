subroutine SetLin_KS(lnr)

use nrtype
use ifc_int_ks

implicit none

real(dp), dimension(:), intent(out) :: lnr
! Returns the linear operator Lin in Fourier space for KSe. 

integer(i4b) :: k

lnr = 0.0_dp

do k=1,size(lnr)
	lnr(k) = ((Real(k-1)/L)**2-(Real(k-1,dp)/L)**4)
end do

end subroutine
