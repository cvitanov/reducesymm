module parameters

use nrtype

implicit none

real(dp), parameter :: Tguess=5.881088923425644, alpha=0.2_dp, beta=0.2_dp, gamma=5.7_dp
real(dp), parameter :: tol=1e-13, Delta_x=0.001_dp, direction=1.0_dp 
integer(i4b), parameter ::  d=3, sect=1, maxIter=1000, nstepsN = 10000, nsteps = 100000, nstepsP=1000000
integer(i4b), parameter :: nInters=1
CHARACTER(len=*), parameter :: format_label3='(3F17.14)',format_label='(4F17.14)'

INTERFACE
	subroutine init_a(a)
	use nrtype
	Real(dp), dimension(:), intent(inout) :: a
	end subroutine
END INTERFACE

contains

	subroutine init_a(a)
	use nrtype
	Real(dp), dimension(:), intent(inout) :: a
	
	a=0.0_dp
	a(sect)=1.0_dp
	end subroutine

end module
