module ifc_newt

use nrtype

integer(i4b) :: newton_condition_met

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
        SUBROUTINE mnewt_elim(ntrial,elim,x,tolx,tolf,usrfun)
                USE nrtype
                IMPLICIT NONE
                INTEGER(I4B), INTENT(IN) :: ntrial,elim
                REAL(dp), INTENT(IN) :: tolx,tolf
                real(dpc), DIMENSION(:), INTENT(INOUT) :: x
                INTERFACE
                        SUBROUTINE usrfun(x,elim,fvec,fjac)
                        USE nrtype
                        IMPLICIT NONE
                        real(dpc), DIMENSION(:), INTENT(IN) :: x
			integer(i4b), DIMENSION(:), INTENT(IN) :: elim
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
	SUBROUTINE mnewtRPOdamp(ntrial,x,tolx,tolf,T,kappa,damp,usrfun)
		USE nrtype
		IMPLICIT NONE
		INTEGER(I4B), INTENT(IN) :: ntrial
		REAL(dp), INTENT(IN) :: tolx,tolf
		real(dpc), DIMENSION(:), INTENT(INOUT) :: x
		real(dp), intent(inout) :: T,kappa
		real(dp), intent(in):: damp
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
