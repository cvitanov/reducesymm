	SUBROUTINE sort_pick_r(arr)
	USE nrtype
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
	INTEGER(I4B) :: i,j,n
	REAL(DP) :: a
	n=size(arr)
	do j=2,n
		a=arr(j)
		do i=j-1,1,-1
			if (arr(i) <= a) exit
			arr(i+1)=arr(i)
		end do
		arr(i+1)=a
	end do
	END SUBROUTINE sort_pick_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE sort_pick_Re(arr)
	USE nrtype
	IMPLICIT NONE
	complex(dpc), DIMENSION(:), INTENT(INOUT) :: arr
	INTEGER(I4B) :: i,j,n
	complex(dpc) :: a
	!! sort_pick_real sorts complex number comparing their real parts
	n=size(arr)
	do j=2,n
		a=arr(j)
		do i=j-1,1,-1
			if (real(arr(i)) <= real(a)) exit
			arr(i+1)=arr(i)
		end do
		arr(i+1)=a
	end do
	END SUBROUTINE sort_pick_Re