!     Last change:  V    17 Mar 2004    1:52 am
function roesslerVar(x,y)

USE nrtype
USE parameters

IMPLICIT NONE

REAL(DP), INTENT(IN) :: x
REAL(DP), DIMENSION(:), INTENT(IN) :: y

REAL(DP), DIMENSION(size(y),size(y)) :: roesslerVar

roesslerVar(1,1)=0
roesslerVar(1,2)=-1
roesslerVar(1,3)=-1
roesslerVar(2,1)=1
roesslerVar(2,2)=alpha
roesslerVar(2,3)=0
roesslerVar(3,1)=y(3)
roesslerVar(3,2)=0
roesslerVar(3,3)=y(1)-gamma



end function roesslerVar