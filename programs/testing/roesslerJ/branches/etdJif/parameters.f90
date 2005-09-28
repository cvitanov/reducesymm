module parameters

use nrtype

implicit none

complex(dpc), parameter :: alpha=(0.2,0.0), beta=(0.2,0.0), gamma=(5.7,0.0)
integer(i4b), parameter ::  d=3 
CHARACTER(len=*), parameter :: format_label='(3F16.10)'

end module