module ifc_newt

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