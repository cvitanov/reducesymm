	SUBROUTINE zbracB(func,x1,x2,expectedTrue,expectedFalse,succes)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(DP), INTENT(INOUT) :: x1,x2
	LOGICAL(LGT), INTENT(OUT) :: succes
	integer(i4b), intent(in) :: expectedTrue,expectedFalse
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		LOGICAL(LGT) :: func
		END FUNCTION func
	END INTERFACE
	! Brackets a point where the outcome of a Boolean function changes.
	INTEGER(I4B), PARAMETER :: NTRY=50
	REAL(DP), PARAMETER :: FACTOR=1.6_DP
	INTEGER(I4B) :: j
	LOGICAL(LGT) :: f(1:2)
	if (x1 == x2) call nrerror('zbrac: you have to guess an initial range')
	f(1)=func(x1)
	f(2)=func(x2)
	succes=.true.
	do j=1,NTRY
		if ( f(expectedTrue) .and. (.not. f(expectedFalse)) ) return
		if ( .not. f(expectedTrue) ) then
			if ( expectedTrue == 1 ) then
				x1=x1+FACTOR*(x1-x2)
				f(1)=func(x1)
			else
				x2=x2+FACTOR*(x2-x1)
				f(2)=func(x2)
			end if
		else
			if ( expectedFalse == 1 ) then
				x1=x1+FACTOR*(x1-x2)
				f(1)=func(x1)
			else
				x2=x2+FACTOR*(x2-x1)
				f(2)=func(x2)
			end if		
		end if
	end do
	succes=.false.
	END SUBROUTINE zbracB
