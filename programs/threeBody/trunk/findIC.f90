subroutine findIC(xi,p_xi,E,Z)

use nrtype

implicit none

real(dp), dimension(:), intent(inout) :: xi
real(dp), dimension(:), intent(in) :: p_xi
real(dp), intent(in) :: E
real(dp), intent(in) :: Z
!!
real(dp), dimension(size(xi)):: xidum
real(dp) :: alpha,beta,gamma

!!! Find Abs(xi(2)) compatible with E
alpha = E-0.5_dp*p_xi(1)**2+Z/xi(1)
beta = (E-0.5_dp*p_xi(1)**2)*xi(1)+2*Z-1
gamma = Z*xi(1)
xi(2) = (-beta-sqrt(beta**2-4*alpha*gamma))/(2*alpha)
xidum(2) = (-beta+sqrt(beta**2-4*alpha*gamma))/(2*alpha)
if (xi(2)<0) xi(2)=xidum(2)
xi(2)=-xi(2) ! xi(2) is negative for us.
!print *,xi(1),xi(2)
! The energy value calculated using the values we have found is slightly different. Take care of this?
!E=0.5_dp*p_xi(1)**2+0.5_dp*p_xi(2)**2-Z/Abs(xi(1))-Z/Abs(xi(2))+1/Sqrt((xi(1)-xi(2))**2)


end subroutine findIC
