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

dum=assert_eq(size(ANdiag,1),size(ANdiag,2), size(a)-1, 'SetANDiag_KS')

ANdiag=(0.0,0.0)

do k=1,size(ANdiag,1)
	do j=k,size(ANdiag,2)
!		if ( j+k .le. size(ANdiag,2) ) then
!			ANdiag(k,j) = ii*k*(conjg(a(Abs(k-j)))+a(k+j))/L
!		else
			ANdiag(k,j) = 2*ii*k*(conjg(a(Abs(k-j))))/L
!		endif
	end do
	do j=1,k-1
!		if ( j+k .le. size(ANdiag,2) ) then
!			ANdiag(k,j) = ii*k*(a(k-j)+a(k+j))/L
!		else
			ANdiag(k,j) = 2*ii*k*(a(k-j))/L
!		endif
	end do
end do

end subroutine
