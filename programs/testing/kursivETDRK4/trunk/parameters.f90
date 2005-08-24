module parameters

use nrtype

real(dp), parameter :: L=6.51739492_dp ! Dimensionless length used here
integer(i4b), parameter :: M=32, d=256
real(dp), parameter :: h=0.025_dp, R=1.0_dp, ti=0.0_dp,tf=2000.0_dp
integer(i4b),parameter :: Nplt=2000
character(len=*), parameter :: frm_t='(1000F12.4)', frm_u='(256F19.12)', frm_a='(129F19.12)'

end module parameters
