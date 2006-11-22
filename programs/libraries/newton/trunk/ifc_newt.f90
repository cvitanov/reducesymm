module ifc_newt

use nrtype

integer(i4b) :: newton_condition_met=0

interface
	SUBROUTINE mnewt(ntrial,x,tolx,tolf,usrfun)
		USE nrtype
		IMPLICIT NONE
		INTEGER(I4B), INTENT(IN) :: ntrial
		REAL(dp), INTENT(IN) :: tolx,tolf
		real(dpc), DIMENSION(:), INTENT(INOUT) :: x
		INTERFACE
			SUBROUTINE usrfun(x,fvec,fjac)
			USE nrtype
			IMPLICIT NONE
			real(dpc), DIMENSION(:), INTENT(IN) :: x
			real(dpc), DIMENSION(:), INTENT(OUT) :: fvec
			real(dpc), DIMENSION(:,:), INTENT(OUT) :: fjac
			END SUBROUTINE usrfun
		END INTERFACE
	end subroutine
end interface 

interface
	SUBROUTINE mnewtPO(ntrial,x,tolx,tolf,T,q,usrfun)
		USE nrtype
		IMPLICIT NONE
		INTEGER(I4B), INTENT(IN) :: ntrial
		REAL(dp), INTENT(IN) :: tolx,tolf
		real(dpc), DIMENSION(:), INTENT(INOUT) :: x
		real(dp), intent(inout) :: T
		real(dp), dimension(:), intent(in):: q
		INTERFACE
			SUBROUTINE usrfun(x,fvec,fjac,T,q)
			USE nrtype
			IMPLICIT NONE
			real(dpc), DIMENSION(:), INTENT(IN) :: x
			real(dpc), DIMENSION(:), INTENT(OUT) :: fvec
			real(dpc), DIMENSION(:,:), INTENT(OUT) :: fjac
			real(dp), intent(in) :: T
			real(dp), dimension(:), intent(in) :: q
			END SUBROUTINE usrfun
		END INTERFACE
	end subroutine
end interface 

interface
	SUBROUTINE mnewtRPO(ntrial,x,tolx,tolf,T,kappa,usrfun)
		USE nrtype
		IMPLICIT NONE
		INTEGER(I4B), INTENT(IN) :: ntrial
		REAL(dp), INTENT(IN) :: tolx,tolf
		real(dpc), DIMENSION(:), INTENT(INOUT) :: x
		real(dp), intent(inout) :: T,kappa
		INTERFACE
			SUBROUTINE usrfun(x,fvec,fjac,T,kappa)
			USE nrtype
			IMPLICIT NONE
			real(dpc), DIMENSION(:), INTENT(IN) :: x
			real(dpc), DIMENSION(:), INTENT(OUT) :: fvec
			real(dpc), DIMENSION(:,:), INTENT(OUT) :: fjac
			real(dp), intent(in) :: T,kappa
			END SUBROUTINE usrfun
		END INTERFACE
	end subroutine
end interface 

interface
	SUBROUTINE mnewtRPOdamp(ntrial,x,tolx,tolf,T,kappa,q,damp,usrfun)
		USE nrtype
		IMPLICIT NONE
		INTEGER(I4B), INTENT(IN) :: ntrial
		REAL(dp), INTENT(IN) :: tolx,tolf
		real(dpc), DIMENSION(:), INTENT(INOUT) :: x
		real(dp), intent(inout) :: T,kappa
		real(dp), dimension(:), intent(in):: q
		real(dp), intent(in):: damp
		INTERFACE
			SUBROUTINE usrfun(x,fvec,fjac,T,kappa,q)
			USE nrtype
			IMPLICIT NONE
			real(dpc), DIMENSION(:), INTENT(IN) :: x
			real(dpc), DIMENSION(:), INTENT(OUT) :: fvec
			real(dpc), DIMENSION(:,:), INTENT(OUT) :: fjac
			real(dp), intent(in) :: T,kappa
			real(dp), dimension(:), intent(in) :: q
			END SUBROUTINE usrfun
		END INTERFACE
	end subroutine
end interface 

interface
	SUBROUTINE mnewtTW(ntrial,x,tolx,tolf,kappa,usrfun)
		USE nrtype
		IMPLICIT NONE
		INTEGER(I4B), INTENT(IN) :: ntrial
		REAL(dp), INTENT(IN) :: tolx,tolf
		real(dpc), DIMENSION(:), INTENT(INOUT) :: x
		real(dp), intent(inout) :: kappa
		INTERFACE
			SUBROUTINE usrfun(x,fvec,fjac,kappa)
			USE nrtype
			IMPLICIT NONE
			real(dpc), DIMENSION(:), INTENT(IN) :: x
			real(dpc), DIMENSION(:), INTENT(OUT) :: fvec
			real(dpc), DIMENSION(:,:), INTENT(OUT) :: fjac
			real(dp), intent(in) :: kappa
			END SUBROUTINE usrfun
		END INTERFACE
	end subroutine
end interface 

interface
	SUBROUTINE mnewt_c(ntrial,x,tolx,tolf,usrfun)
		USE nrtype
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
	end subroutine
end interface 

end module