	SUBROUTINE mnewt_c(ntrial,x,tolx,tolf,usrfun)
	USE nrtype
	use la_precision, only: wp => dp
	use f95_lapack, only: LA_GESV
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ntrial
	REAL(dp), INTENT(IN) :: tolx,tolf
	complex(dpc), DIMENSION(:), INTENT(INOUT) :: x
	INTERFACE
		SUBROUTINE usrfun(x,fvec,fjac)
		USE nrtype
		IMPLICIT NONE
		complex(dpc), DIMENSION(:), INTENT(IN) :: x
		complex(dpc), DIMENSION(:), INTENT(OUT) :: fvec
		complex(dpc), DIMENSION(:,:), INTENT(OUT) :: fjac
		END SUBROUTINE usrfun
	END INTERFACE
	INTEGER(I4B) :: i
	complex(dpc), DIMENSION(size(x)) :: fvec,p
	complex(dpc), DIMENSION(size(x),size(x)) :: fjac
	do  i=1,ntrial
		call usrfun(x,fvec,fjac)
		if (sum(abs(fvec)) <= tolf) RETURN
		p=-fvec
		call la_gesv(fjac,p)
		x=x+p
		if (sum(abs(p)) <= tolx) RETURN
	end do
	END SUBROUTINE mnewt_c
