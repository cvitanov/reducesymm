subroutine KSMcGtoCart_TBL(Et,Q,P,hyperR,x,p_x,E)

use nrtype

implicit none

real(dp), dimension(:), intent(in) :: Q,P
real(dp), intent(in)::  Et, hyperR
real(dp), dimension(:), intent(out):: x,p_x
real(dp), intent(out):: E 

real(dp), dimension(size(Q)) :: Qdum,rdum

!
x(1)=hyperR*Q(1)**2
x(2)=-hyperR*Q(2)**2
p_x(1)=Q(1)*P(1)/(2*Sqrt(hyperR)*Q(1)**2)
p_x(2)=-Q(2)*P(2)/(2*Sqrt(hyperR)*Q(2)**2)
E=Et/hyperR

end subroutine