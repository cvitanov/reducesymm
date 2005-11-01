subroutine SetANdiag_KS(a,ANdiag)
use nrtype
use parameters
implicit none
complex(dpc), dimension(:), intent(in) :: a
complex(dpc), dimension(:,:), intent(out) :: ANdiag
!
!
integer(i4b) :: k,j

do k=1,size(a)
	do j=1,size(ANdiag,1)
		ANdiag(k,j) = 2*ii*(k-1)*a(1+Abs(k-j))/L
	end do
end do

end subroutine
