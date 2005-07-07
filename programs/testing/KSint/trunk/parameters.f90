module parameters

use nrtype

implicit none

! L is the dimensionless length
real(dp), parameter :: nudum=0.021_dp !L = 6.900655593_dp 
real(dp), parameter :: L=6.9006555934235427330_dp
CHARACTER(len=*), parameter :: format_label='(32F12.8)', flu='(150F12.8)'
complex(dpc), parameter :: ii=(0.0_dp,1.0_dp)

end module
