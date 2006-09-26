MODULE ifc_util

INTERFACE
	FUNCTION UnitMatrix(N)
		USE nrtype
		IMPLICIT NONE
		INTEGER(I4B) :: N
		REAL(DP) :: UnitMatrix(N,N)
	END FUNCTION UnitMatrix
END INTERFACE


interface DiagMul
	function DiagMul_rc(v,M,N)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(in) :: v
		complex(dpc), dimension(:,:), intent(in) :: M
		integer(i4b), intent(in) :: N
		complex(dpc), dimension(N,N) :: DiagMul_rc
	end function DiagMul_rc
	function DiagMul_cc(v,M,N)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in) :: v
		complex(dpc), dimension(:,:), intent(in) :: M
		integer(i4b), intent(in) :: N
		complex(dpc), dimension(N,N) :: DiagMul_cc
	end function DiagMul_cc
	function DiagMul_rr(v,M,N)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(in) :: v
		real(dp), dimension(:,:), intent(in) :: M
		integer(i4b), intent(in) :: N
		real(dp), dimension(N,N) :: DiagMul_rr
	end function DiagMul_rr
end interface 

interface
	subroutine donothing(tmp)
		use nrtype
		implicit none
		real(dpc), intent(inout):: tmp
	end subroutine
end interface

interface
	logical function SelectLargeEig(wR_j,wI_j)
		use nrtype
		implicit none
		real(dp), INTENT(IN) :: wR_j,wI_j
	end function
end interface

interface
	logical function SelectSmallEig(wR_j,wI_j)
		use nrtype
		implicit none
		real(dp), INTENT(IN) :: wR_j,wI_j
	end function
end interface

interface sort_pick
	SUBROUTINE sort_pick_r(arr)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
	END SUBROUTINE sort_pick_r
	SUBROUTINE sort_pick_Re(arr)
		USE nrtype
		IMPLICIT NONE
		complex(dpc), DIMENSION(:), INTENT(INOUT) :: arr
	END SUBROUTINE sort_pick_Re
end interface

END MODULE