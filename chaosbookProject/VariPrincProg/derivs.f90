!     Last change:  V    20 Feb 2004    2:24 am
SUBROUTINE derivs(t,y,dydt)
USE nrtype
IMPLICIT NONE
REAL(DP), INTENT(IN) :: t
REAL(DP), DIMENSION(:), INTENT(IN) :: y
REAL(DP), DIMENSION(:), INTENT(OUT) :: dydt

! Derivs function for KS system. NOTICE: \nu is input here.



REAL(DP) :: nu
REAL(DP) :: sm
INTEGER :: d
INTEGER :: m,k




 nu=0.029910_DP

 d = size(y) 

  DO k=1,d

   sm=0_dp

   DO m=1,k
     if (k-m .NE. 0) sm=sm+y(m)* Sign(1,k-m) * y(Abs(k-m))
   END DO

   DO m=k+1,d
     sm=sm-y(m)* Sign(1,m-k) * y(Abs(m-k)) 
   END DO

   DO m=1,d-k
	 sm=sm-y(m)*y(k+m)
   END DO

   dydt(k)= y(k)*(1-nu*k**2)*k**2 - k*sm 

 END DO


END SUBROUTINE derivs
