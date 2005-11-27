subroutine SetLin_KS(Lin)

use nrtype
use parameters

implicit none

real(dp), dimension(:), intent(out) :: Lin
! Returns the linear operator Lin in Fourier space for KSe. 

integer(i4b) :: k

Lin = 0.0_dp

do k=1,size(Lin)
	Lin(k) = (((k-1)/L)**2-((k-1)/L)**4)
end do

end subroutine
