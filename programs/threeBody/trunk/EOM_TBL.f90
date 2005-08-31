subroutine EOM_TBL(tau,f,v)

use nrtype
use parameters

implicit none

real(dp), dimension(:), intent(in) :: f
real(dp), intent(in) :: tau
real(dp), dimension(:), intent(out) :: v
 
! All quantities are post McGehee. 
! Change of vector indices with respect to the references: 4 --> 2
! Yet r_1, r_2 agree with references.

real(dp), dimension(2) :: Q, P, r
real(dp) :: E, hyperR
real(dp) :: p_R,r_12

! Make code readable:
Q=f(1:2)
P=f(3:4)
E=f(5)
!print *,"E?",E
hyperR=f(6)
r=Q**2

p_R=(Q(1)*P(1)+Q(2)*P(2))/2
r_12=Sqrt(Q(1)**4+Q(2)**4+2*(Q(1)*Q(2))**2)

!print *,"r_12", r_12

!dQ1/dtau 
v(1)= r(2)*P(1)/4.0_dp-Q(1)*r(1)*r(2)*p_R/2.0_dp
!dQ4/dtau
v(2)= r(1)*P(2)/4.0_dp-Q(2)*r(1)*r(2)*p_R/2.0_dp
!dP1/dtau
v(3)= -( (1/4)*Q(1)*p(2)**2 - 2*Z*Q(1) + 2*Q(1)*r(2)*(-E+1/(r_12))-2*r(1)*r(2)*(r(1)*Q(1)+Q(1)*Q(2)**2)/r_12**3  )
!dP4/dtau
v(4)= -( (1/4)*Q(2)*p(1)**2 - 2*Z*Q(2) + 2*Q(2)*r(1)*(-E+1/(r_12))-2*r(1)*r(2)*(r(2)*Q(2)+Q(2)*Q(1)**2)/r_12**3  )
!dE/dtau
v(5)= r(1)*r(2)*p_R*E
!print *,"dE/dt",v(5)
!dR/dtau
v(6)= hyperR*r(1)*r(2)*p_R
!dt/dtau 
v(7)= r(1)*r(2)*Sqrt(hyperR)**3 

end subroutine