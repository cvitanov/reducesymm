subroutine SetANdiag_KS(a,ANdiag)
use nrtype
use nrutil
use parameters
implicit none
complex(dpc), dimension(:), intent(in) :: a
complex(dpc), dimension(:,:), intent(out) :: ANdiag
!
!
integer(i4b) :: k,j,dum

dum=assert_eq(size(ANdiag,1),size(ANdiag,2), size(a), 'SetANDiag_KS')

ANdiag=(0.0,0.0)

do k=0,size(ANdiag,1)-1
	do j=k,size(ANdiag,2)-1
		if ( j+k .le. size(ANdiag,2)-1 ) then
			ANdiag(k+1,j+1) = ii*k*(conjg(a(1+Abs(k-j)))+a(1+k+j))/L
		else
			ANdiag(k+1,j+1) = ii*k*(conjg(a(1+Abs(k-j))))/L
		endif
	end do
	do j=0,k-1
		if ( j+k .le. size(ANdiag,2)-1 ) then
			ANdiag(k+1,j+1) = ii*k*(a(1+k-j)+a(1+k+j))/L
		else
			ANdiag(k+1,j+1) = ii*k*(a(1+k-j))/L
		endif
	end do
end do

end subroutine
