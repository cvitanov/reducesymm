function Rc_M(w,dmn)

use nrtype
use nrutil, only: assert_eq

implicit none

integer(i4b), intent(in) :: dmn
complex(dpc), dimension(dmn) :: Rc_M
real(dp), intent(in) :: w
!! Returns the matrix the action of which performs rotation in 
!! each imaginary plane by angle w. This
!! corresponds to a translation by kappa=w*Ltilde. 
integer(i4b):: k,dum

dum=assert_eq(dmn,size(Rc_M),"Rc_M")

Rc_M=(0.0_dp,0.0_dp)

do k=1,dmn
	Rc_M(k)=cdExp(ii*Real(k-1,dp)*w)
end do

end function Rc_M


function Rc_a(w,a,dmn)

use nrtype
use nrutil, only: assert_eq

implicit none

integer(i4b), intent(in) :: dmn
complex(dpc), dimension(dmn) :: Rc_a
complex(dpc), dimension(dmn), intent(in) :: a 
real(dpc), intent(in) :: w
!! Performs rotation in each imaginary plane by angle w. This
!! corresponds to a translation by kappa=w*Ltilde. 
integer(i4b):: k,dum

dum=assert_eq(dmn,size(Rc_a),size(a),"Rc_a")

Rc_a=(0.0_dp,0.0_dp)

do k=1,dmn
	Rc_a(k)=cdExp(ii*Real(k-1,dp)*w)*a(k)
end do

end function Rc_a