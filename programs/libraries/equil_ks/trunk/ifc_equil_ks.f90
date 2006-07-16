module ifc_equil_ks

use nrtype

integer(i4b) :: ntrial
integer(i4b) :: d 
real(dp) :: L
real(dp) :: tolbc, tolf

!!!!
interface
	subroutine ksFJ(bc,fvec,fjac)
		USE nrtype
		implicit none
		real(dp), DIMENSION(:), INTENT(IN) :: bc
		real(dp), DIMENSION(:), INTENT(OUT) :: fvec
		real(dp), DIMENSION(:,:), INTENT(OUT) :: fjac
	end subroutine
end interface
!!!!

interface
	subroutine SetNlin_KS(a,N_a)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in) :: a
		complex(dpc), dimension(:), intent(out) :: N_a
	end subroutine
end interface

end module
