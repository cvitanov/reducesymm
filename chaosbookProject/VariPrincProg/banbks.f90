SUBROUTINE banbks(a,m1,m2,al,indx,b)

! From NUMERICAL RECIPES IN FORTRAN 90: THE Art of PARALLEL Scientific Computing

USE nrtype; USE nrutil, ONLY : assert_eq,swap
IMPLICIT NONE

REAL(DP), DIMENSION(:,:), INTENT(IN) :: a,al
INTEGER(I4B), INTENT(IN) :: m1,m2
INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
!Given the arrays a, al, and indx as returned from bandec, and given a right-hand-side
!vector b, solves the band diagonal linear equations A·x = b. The solution vector x overwrites
!b. The other input arrays are not modified, and can be left in place for successive calls with
!different right-hand sides.

INTEGER(I4B) :: i,k,l,mdum,mm,n

n=assert_eq(size(a,1),size(al,1),size(b),size(indx),'banbks: n')
mm=assert_eq(size(a,2),m1+m2+1,'banbks: mm')
mdum=assert_eq(size(al,2),m1,'banbks: mdum')

do k=1,n !Forward substitution, unscrambling the permuted rows as we go. 
	l=min(n,m1+k)
	i=indx(k)
	if (i /= k) call swap(b(i),b(k))
	b(k+1:l)=b(k+1:l)-al(k,1:l-k)*b(k)
end do
do i=n,1,-1 !Backsubstitution.
	l=min(mm,n-i+1)
	b(i)=(b(i)-dot_product(a(i,2:l),b(1+i:i+l-1)))/a(i,1)
end do

END SUBROUTINE banbks