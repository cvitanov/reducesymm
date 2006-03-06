module ifc_newton

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