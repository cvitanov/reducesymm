subroutine TBL_CartToKSMcG(x,p_x,dt,E,Q,P,r,dtau,Et,R)

use nrtype

implicit none

real(dp), dimension(:), intent(in):: x,p_x
real(dp), intent(in):: dt, E 
real(dp), dimension(:), intent(out) :: Q,P,r
real(dp), intent(out):: dtau, Et, R


!Part I: Kustaanheimo-Stiefel
Q(1)=Sqrt(x(1))
Q(2)=-Sqrt(x(2))
hyperR=Sqrt(Q(1)**4+Q(2)**4)
r(1)=Q(1)**2
r(2)=Q(2)**2
P(1)=2*r(1)*p_x(1)/Q(1)
P(2)=-2*r(2)*p_x(2)/Q(2)

dtau=dt/(r(1)*r(2))

!Part II: McGehee

Q=Q/Sqrt(hyperR)
dtau=Sqrt(hyperR)*dtau
r=r/hyperR
Et=hyperR*E


end subroutine