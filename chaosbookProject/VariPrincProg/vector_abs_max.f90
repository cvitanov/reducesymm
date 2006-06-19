FUNCTION vector_abs_max(b)

USE nrtype
IMPLICIT NONE

REAL(DP) :: vector_abs_max
REAL(DP), INTENT(IN) :: b(:)

!Returns the maximum absolute value of the elements of a vector b.


INTEGER(I4B) i


vector_abs_max=Abs(b(1))

do i=1,size(b)
	 if (Abs(b(i)) > vector_abs_max) vector_abs_max=Abs(b(i))
end do


END FUNCTION