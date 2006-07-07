	SUBROUTINE mnewtRPO(ntrial,x,tolx,tolf,T,kappa,usrfun)
	USE nrtype
	use la_precision, only: wp => dp
	use f95_lapack, only: LA_GESV
	use nrutil, only: assert_eq
	use ifc_newt
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ntrial
	REAL(DP), INTENT(IN) :: tolx,tolf
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
	real(dp), intent(inout) :: T,kappa
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

	do  i=1,ntrial
		print *,"Newton iteration #", i
		call usrfun(x,fvec,fjac,T,kappa)
		if (sum(abs(fvec)) <= tolf) then
			print *,"Condition for fvec met, with fvec=",sum(abs(fvec))
			newton_condition_met=1
			RETURN
		endif
		p=-fvec
		call la_gesv(fjac,p)
		x=x+p(1:size(p)-2)
		T=T+p(size(p)-1)
		kappa=kappa+p(size(p))
		print *,"DT",p(size(p)-1),p(size(p))
		print *,"Delta p",sum(abs(p))
		if (sum(abs(p)) <= tolx) then 
			print *,"Condition for tolx met, with tolx=",sum(abs(p))
			newton_condition_met=2
			RETURN
		end if
	end do

	newton_condition_met=0

	END SUBROUTINE mnewtRPO
