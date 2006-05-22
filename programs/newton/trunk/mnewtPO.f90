	SUBROUTINE mnewtPO(ntrial,x,tolx,tolf,T,q,usrfun)
	USE nrtype
	use la_precision, only: wp => dp
	use f95_lapack, only: LA_GESV
	use nrutil, only: assert_eq
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ntrial
	REAL(DP), INTENT(IN) :: tolx,tolf
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
	real(dp), intent(inout) :: T
	real(dp), dimension(:), intent(in):: q
	INTERFACE
		SUBROUTINE usrfun(x,fvec,fjac,T,q)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(OUT) :: fvec
		REAL(DP), DIMENSION(:,:), INTENT(OUT) :: fjac
		real(dp), intent(in) :: T
		real(dp), dimension(:), intent(in):: q
		END SUBROUTINE usrfun
	END INTERFACE
	INTEGER(I4B) :: i
	REAL(DP) :: ndum
	REAL(DP), DIMENSION(size(x)+1) :: fvec,p
	REAL(DP), DIMENSION(size(x)+1,size(x)+1) :: fjac

	ndum=assert_eq(size(x),size(q),'mnewtPO')

	do  i=1,ntrial
		print *,i
		call usrfun(x,fvec,fjac,T,q)
		if (sum(abs(fvec)) <= tolf) then
			RETURN
		endif
		p=-fvec
		call la_gesv(fjac,p)
		x=x+p(1:size(p)-1)
		T=T+p(size(p))
		if (sum(abs(p)) <= tolx) then 
			RETURN
		end if
	end do
	END SUBROUTINE mnewtPO
