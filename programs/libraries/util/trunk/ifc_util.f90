MODULE ifc_util

INTERFACE
	FUNCTION UnitMatrix(N)
		USE nrtype
		IMPLICIT NONE
		INTEGER(I4B) :: N
		REAL(DP) :: UnitMatrix(N,N)
	END FUNCTION
END INTERFACE


interface DiagMul
	function DiagMul_rc(v,M,N)
		use nrtype
		implicit none
		real(dp), dimension(N), intent(in) :: v
		complex(dpc), dimension(N,N), intent(in) :: M
		integer(i4b), intent(in) :: N
		complex(dpc), dimension(N,N) :: DiagMul_rc
	end function DiagMul_rc
	function DiagMul_cc(v,M,N)
		use nrtype
		implicit none
		complex(dpc), dimension(N), intent(in) :: v
		complex(dpc), dimension(N,N), intent(in) :: M
		integer(i4b), intent(in) :: N
		complex(dpc), dimension(N,N) :: DiagMul_cc
	end function DiagMul_cc
end interface 



END MODULE