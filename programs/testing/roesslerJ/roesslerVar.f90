!     Last change:  V    17 Mar 2004    1:52 am
function roesslerVar(x,y)

USE nrtype
USE parameters

IMPLICIT NONE

REAL(DP), INTENT(IN) :: x
REAL(DP), DIMENSION(:), INTENT(IN) :: y

REAL(DP), DIMENSION(size(y),size(y)) :: MatVar

MatVar(1,1)=0
MatVar(1,2)=-1
MatVar(1,3)=-1
MatVar(2,1)=1
MatVar(2,2)=alpha
MatVar(2,3)=0
MatVar(3,1)=y(3)
MatVar(3,2)=0
MatVar(3,3)=y(1)-gamma



end function roesslerVar