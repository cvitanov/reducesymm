module parameters

use nrtype

real(dp), parameter :: L=3.0_dp ! Dimensionless length used here
integer(i4b), parameter :: M=16, d=128, maxIter=1000
real(dp), parameter :: h=0.025_dp, R=1.0_dp
real(dp), parameter :: tol=1e-8
character(len=*), parameter :: frm_t='(1000F12.4)', frm_u='(128F19.12)'

end module parameters
