subroutine vectorize(Q,P,Et,hyperR,ti,fi)

use nrtype
use nrutil, only: assert_eq

implicit none

real(dp), dimension(:), intent(in):: Q,P
real(dp), intent(in):: Et,hyperR,ti
real(dp), dimension(:), intent(out):: fi 
!
integer :: ndum

ndum = assert_eq(size(fi),7,'vectorize')

fi(1:2)=Q
fi(3:4)=P
fi(5)= Et
fi(6)=hyperR
fi(7)=ti


end subroutine vectorize

