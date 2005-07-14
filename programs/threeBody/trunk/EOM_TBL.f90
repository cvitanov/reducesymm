subroutine EOM_TBL(Q,P,r,E,hyperR,v)

use nrtype
use parameters

implicit none

real(dp), dimension(:), intent(in) :: Q,P,r
real(dp), intent(in) :: E,hyperR
real(dp), dimension(:), intent(out) :: v
 
! All quantities are post McGehee. 
! Change of vector indices with respect to the references: 4 --> 2
! Yet r_1, r_2 agree with references.

real(dp) :: p_R,r_12

p_R=(Q(1)*P(1)+Q(2)*P(2))/2
r_12=Abs(Q(1)+Q(2))

!dQ1/dtau 
v(1)= r(2)*P(1)/4.0_dp-Q(1)*r(1)*r(2)*p_R/2.0_dp
!dQ4/dtau
v(2)= r(1)*P(2)/4.0_dp-Q(2)*r(1)*r(2)*p_R/2.0_dp
!dP1/dtau
v(3)= -( (Q(1)*p(2)**2)/4 - 2*Z*Q(1) + 2*Q(1)*r(2)*(-E+1/(r_12))-2*r(1)*r(2)*(r(1)*Q(1)+Q(2)**2*Q(1))/r_12**3  )
!dP4/dtau
v(4)= -( (Q(2)*p(1)**2)/4 - 2*Z*Q(2) + 2*Q(2)*r(1)*(-E+1/(r_12))-2*r(1)*r(2)*(r(2)*Q(2)+Q(1)**2*Q(2))/r_12**3  )
!dE/dtau
v(5)= r(1)*r(2)*p_R*E
!dR/dtau
v(6)= hyperR*r(1)*r(2)*p_R
!dt/dtau 
v(7)= r(1)*r(2)*Sqrt(hyperR)**3 

end subroutine