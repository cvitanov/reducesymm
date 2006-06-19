SUBROUTINE banmul(a,m1,m2,x,b)

! From NUMERICAL RECIPES IN FORTRAN 90: THE Art of PARALLEL Scientific Computing

USE nrtype; USE nrutil, ONLY : assert_eq,arth
IMPLICIT NONE
REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
INTEGER(I4B), INTENT(IN) :: m1,m2
REAL(DP), DIMENSION(:), INTENT(IN) :: x
REAL(DP), DIMENSION(:), INTENT(OUT) :: b
!Matrix multiply b = A · x, where A is band diagonal with m1 rows below the diagonal and
!m2 rows above. If the input vector x and output vector b are of length N, then the array
!a(1:N,1:m1+m2+1) stores A as follows: The diagonal elements are in a(1:N,m1+1).
!Subdiagonal elements are in a(j:N,1:m1) (with j > 1 appropriate to the number of
!elements on each subdiagonal). Superdiagonal elements are in a(1:j,m1+2:m1+m2+1)
!with j < N appropriate to the number of elements on each superdiagonal.
INTEGER(I4B) :: m,n
REAL(DP), DIMENSION(size(a,1),size(a,2)) :: y
n=assert_eq(size(a,1),size(b),size(x),'banmul: n')
m=assert_eq(size(a,2),m1+m2+1,'banmul: m')
y=spread(x,dim=2,ncopies=m) !Duplicate x into columns of y
y=eoshift(y,dim=1,shift=arth(-m1,1,m)) !Shift columns by a linear progression
b=sum(a*y,dim=2) !Multiply by the band-diagonal elements, and sum.

END SUBROUTINE banmul