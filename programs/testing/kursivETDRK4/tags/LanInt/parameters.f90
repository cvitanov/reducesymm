module parameters

use nrtype

real(dp), parameter :: L=6.51739491961311490087_dp ! Dimensionless length used here
integer(i4b), parameter :: M=32, d=32, Niter=5
real(dp), parameter :: h=0.0001_dp, R=1.0_dp, ti=0.0_dp,tf=100.0_dp
integer(i4b),parameter :: Nplt=1000
character(len=*), parameter :: frm_t='(100F12.4)', frm_u='(32F19.12)', frm_a='(17F19.12)'

end module parameters
