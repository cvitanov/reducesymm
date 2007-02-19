	SUBROUTINE mnewtRPOdamp(ntrial,x,tolx,tolf,T,kappa,damp,usrfun)
	USE nrtype
	use la_precision, only: wp => dp
	use f95_lapack, only: LA_GESV
	use nrutil, only: assert_eq
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ntrial
	REAL(DP), INTENT(IN) :: tolx,tolf
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
	real(dp), intent(inout) :: T,kappa
	real(dp), intent(in):: damp
	INTERFACE
		SUBROUTINE usrfun(x,fvec,fjac,T,kappa)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(OUT) :: fvec
		REAL(DP), DIMENSION(:,:), INTENT(OUT) :: fjac
		real(dp), intent(in) :: T, kappa
		END SUBROUTINE usrfun
	END INTERFACE
	INTEGER(I4B) :: i
	REAL(DP) :: ndum
	REAL(DP), DIMENSION(size(x)+2) :: fvec,p
	REAL(DP), DIMENSION(size(x)+2,size(x)+2) :: fjac
	real(dp) :: alpha ,err

	do  i=1,ntrial
		print *,"Newton iteration #", i
		call usrfun(x,fvec,fjac,T,kappa)
		err=sum(abs(fvec))
		if (err <= tolf) then
			RETURN
		endif
		p=-fvec
		call la_gesv(fjac,p)
		alpha=Exp(-damp*err)
		x=x+alpha*p(1:size(p)-2)
		T=T+alpha*p(size(p)-1)
		kappa=kappa+alpha*p(size(p))
		if (sum(abs(p)) <= tolx) then 
			RETURN
		end if
	end do
	END SUBROUTINE mnewtRPOdamp
