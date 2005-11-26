module parameters

use nrtype

real(dp), parameter :: L=3.063732654518985213439735693086163337284_dp ! Dimensionless length used here
integer(i4b), parameter :: M=32, d=128
integer(i4b),dimension(2), parameter :: seed=(/1727,1234/)
real(dp), parameter :: h=0.025_dp, R=1.0_dp, ti=0.0_dp,tf=10000.0_dp
integer(i4b), parameter :: Nplt=10000
character(len=*), parameter :: frm_t='(1000F12.4)', frm_u='(128F19.12)'

end module parameters
