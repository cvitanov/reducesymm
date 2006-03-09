module parameters

use nrtype

real(dp), parameter :: L=6.1274653090379711529_dp ! Dimensionless length used here.
integer(i4b), parameter :: M=32, d=128
integer(i4b),dimension(2), parameter :: seed=(/1727,1234/)
real(dp), parameter :: h=0.01_dp, R=1.0_dp, ti=3222.0_dp,tf=3245_dp
integer(i4b), parameter :: Nplt=2400
character(len=*), parameter :: frm_t='(1000F12.4)', frm_u='(128F19.12)'

end module parameters
