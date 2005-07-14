subroutine KSMcGtoCart_TBL(Et,Q,P,hyperR,x,p_x,E)

use nrtype

real(dp), dimension(:), intent(in) :: Q,P,r
real(dp), intent(in)::  Et, hyperR
real(dp), dimension(:), intent(out):: x,p_x
real(dp), intent(out):: E 

implicit none

x(1)=Q(1)**2
x(2)=-Q(2)**2
p_x(1)=Q(1)*P(1)
p_x(2)=-Q(4)*P(4)
E=Et/hyperR

end subroutine