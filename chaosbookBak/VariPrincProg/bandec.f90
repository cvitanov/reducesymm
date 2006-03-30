SUBROUTINE bandec(a,m1,m2,al,indx,d)

! From NUMERICAL RECIPES IN FORTRAN 90: THE Art of PARALLEL Scientific Computing

USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,swap,arth
IMPLICIT NONE
REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
INTEGER(I4B), INTENT(IN) :: m1,m2
REAL(DP), DIMENSION(:,:), INTENT(OUT) :: al
INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
REAL(DP), INTENT(OUT) :: d
REAL(DP), PARAMETER :: TINY=1.0e-20_DP

!Given an N × N band diagonal matrix A with m1 subdiagonal rows and m2 superdiagonal
!rows, compactly stored in the array a(1:N,1:m1+m2+1) as described below, this routine 
!constructs an LU decomposition of a rowwise permutation of
!A. The diagonal elements are in a(1:N,m1+1). Subdiagonal elements are in a(j:N,1:m1) 
!(with j > 1 appropriate to the number of elements on each subdiagonal). Superdiagonal elements 
!are in a(1:j,m1+2:m1+m2+1) with j < N appropriate to the number of elements on each 
!superdiagonal. The upper triangular matrix replaces a, while the lower triangular matrix is 
!returned in al(1:N,1:m1). indx is an output vector of length N that records the row permutation
!effected by the partial pivoting; d is output as ±1 depending on whether the number of
!row interchanges was even or odd, respectively. This routine is used in combination with
!banbks to solve band-diagonal sets of equations.

INTEGER(I4B) :: i,k,l,mdum,mm,n
REAL(DP) :: dum

n=assert_eq(size(a,1),size(al,1),size(indx),'bandec: n')
mm=assert_eq(size(a,2),m1+m2+1,'bandec: mm')
mdum=assert_eq(size(al,2),m1,'bandec: mdum')
a(1:m1,:)=eoshift(a(1:m1,:),dim=2,shift=arth(m1,-1,m1)) !Rearrange the storage a bit. 
d=1.0
do k=1,n !For each row...
	l=min(m1+k,n)
	i=imaxloc(abs(a(k:l,1)))+k-1 !Find the pivot element.
	dum=a(i,1)
	if (dum == 0.0) a(k,1)=TINY
	!Matrix is algorithmically singular, but proceed anyway with TINY pivot (desirable in some applications).
	indx(k)=i
	if (i /= k) then !Interchange rows.
		d=-d
		call swap(a(k,1:mm),a(i,1:mm))
	end if
	do i=k+1,l !Do the elimination.
		dum=a(i,1)/a(k,1)
		al(k,i-k)=dum
		a(i,1:mm-1)=a(i,2:mm)-dum*a(k,2:mm)
		a(i,mm)=0.0
	end do
end do

END SUBROUTINE bandec