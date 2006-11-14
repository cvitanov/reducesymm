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
	!! sort_pick_real sorts array arr of complex number comparing their real parts.	
	INTEGER(I4B) :: i,j,n
	complex(dpc) :: a

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE sort_pick_2Re(arr,brr)
	USE nrtype
	USE nrutil, only:assert_eq
	IMPLICIT NONE
	complex(dpc), DIMENSION(:), INTENT(INOUT) :: arr
	complex(dpc), dimension(:,:), intent(inout) :: brr 
	!! sort_pick2_real sorts array arr of complex number comparing their real parts.
	!! It also rearranges the columns of matrix brr to the same order.
	INTEGER(I4B) :: i,j,n
	complex(dpc) :: a
	complex(dpc), dimension(size(brr,2)):: b

	n=assert_eq(size(arr),size(brr,1),size(brr,2),"sort_pick2_Re")
	do j=2,n
		a=arr(j)
		b=brr(:,j)
		do i=j-1,1,-1
			if (real(arr(i)) <= real(a)) exit
			arr(i+1)=arr(i)
			brr(:,i+1)=brr(:,i)
		end do
		arr(i+1)=a
		brr(:,i+1)=b
	end do
	END SUBROUTINE sort_pick_2Re