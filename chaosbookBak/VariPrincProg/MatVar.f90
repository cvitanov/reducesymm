FUNCTION MatVar(a,d)

USE nrtype
IMPLICIT NONE

INTEGER(I4B), INTENT(IN) :: d
REAL(DP) :: MatVar(1:d,1:d)
REAL(DP), INTENT(IN) :: a(:)

!  Computes the matrix of variations for the KS system. NOTICE: \nu is input in here.		
!

REAL(DP) :: nu, aminus, aplus
REAL(DP) :: sm
INTEGER :: m,k,i

aminus=0
aplus=0

nu=0.029910_DP

MatVar=0

DO k=1,d
	DO i=1,d
	 IF (k==i) MatVar(k,i)=(1-nu*k**2)*k**2
	 IF ( (ABS(k-i) .LE. d) .AND. (ABS(k-i) .NE. 0) ) aminus=SIGN(1,k-i)*a(ABS(k-i))
	 IF (k+i .LE. d) aplus=a(k+i)
	 MatVar(k,i)=MatVar(k,i)-2*k*(aminus-aplus)
	END DO
END DO


END FUNCTION
