	SUBROUTINE mnewt(ntrial,x,tolx,tolf,usrfun)
	USE nrtype
	use la_precision, only: wp => dp
	use f95_lapack, only: LA_GESV
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ntrial
	REAL(DP), INTENT(IN) :: tolx,tolf
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
	INTERFACE
		SUBROUTINE usrfun(x,fvec,fjac)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(OUT) :: fvec
		REAL(DP), DIMENSION(:,:), INTENT(OUT) :: fjac
		END SUBROUTINE usrfun
	END INTERFACE
	INTEGER(I4B) :: i
	REAL(DP) :: d
	REAL(DP), DIMENSION(size(x)) :: fvec,p
	REAL(DP), DIMENSION(size(x),size(x)) :: fjac
	real(dp) :: sumx,sumf
	do  i=1,ntrial
		call usrfun(x,fvec,fjac)
		sumf=sum(abs(fvec))
		if ( sumf <= tolf) then
			print *,"Condition for fvec met, with fvec=",sumf
			newton_condition_met=1
			RETURN
		endif
		p=-fvec
		call la_gesv(fjac,p)
		x=x+p
		sumx=sum(abs(p))
		print *,"Delta x",sumx
		if (sum(abs(p)) <= tolx) then 
			print *,"Condition for tolx met, with tolx=",sumx
			newton_condition_met=2
			RETURN
		end if
	end do
	END SUBROUTINE mnewt
