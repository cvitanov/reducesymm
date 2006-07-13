function Rr(w,dmn)

use nrtype
use nrutil, only: assert_eq

implicit none

integer(i4b), intent(in) :: dmn
real(dp), dimension(dmn,dmn) :: Rr
real(dp), intent(in) :: w
!
!
integer(i4b) :: k,dum

dum=assert_eq(dmn,size(Rr,1),size(Rr,1),"Rr")

Rr=0.0_dp

do k=1,dmn/2
	Rr(k,k)= dCos(real(k,dp)*w)
end do
Rr(dmn/2+1:dmn,dmn/2+1:dmn)= Rr(1:dmn/2,1:dmn/2)
do k=1,dmn/2
	Rr(dmn/2+k,k) = -dSin(real(k,dp)*w)
end do
Rr(1:dmn/2,dmn/2+1:dmn)=-Rr(dmn/2+1:dmn,1:dmn/2)

end function Rr