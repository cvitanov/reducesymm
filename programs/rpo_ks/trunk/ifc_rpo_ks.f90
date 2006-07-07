module ifc_rpo_ks

use nrtype

real(dp), dimension(:,:), allocatable :: Jac
real(dp), dimension(:), allocatable :: lin,f0,f1,f2,f3,e,e2
real(dp), dimension(:), allocatable :: f0dum,f1dum,f2dum,f3dum,edum,e2dum
real(dp) :: R, L
integer(i4b) :: M, Nsteps, Nplt, Ntrial

interface
	SUBROUTINE ksFJ(bc,fvec,fjac,T,kappa)
		USE nrtype
		IMPLICIT NONE
		real(dp), DIMENSION(:), INTENT(IN) :: bc
		real(dp), DIMENSION(:), INTENT(OUT) :: fvec
		real(dp), DIMENSION(:,:), INTENT(OUT) :: fjac
		real(dp), intent(in) :: T,kappa
	end subroutine
end interface

interface
	SUBROUTINE ksFJ_diff(bc,fvec,fjac,T,kappa)
		USE nrtype
		IMPLICIT NONE
		real(dp), DIMENSION(:), INTENT(IN) :: bc
		real(dp), DIMENSION(:), INTENT(OUT) :: fvec
		real(dp), DIMENSION(:,:), INTENT(OUT) :: fjac
		real(dp), intent(in) :: T,kappa
	end subroutine
end interface


interface
	subroutine SetLin_KS(lnr)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(out) :: lnr
	end subroutine
end interface
interface
	subroutine SetNlin_KS(a,N_a)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in) :: a
		complex(dpc), dimension(:), intent(out) :: N_a
	end subroutine
end interface


INTERFACE
	SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE derivs	
END INTERFACE
interface
	function MatVar(x,y)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(size(y),size(y)) :: MatVar
	end function MatVar
end interface

interface Rc
	function Rc_M(w,dmn)
		use nrtype
		implicit none
		integer(i4b), intent(in) :: dmn
		complex(dpc), dimension(dmn) :: Rc_M
		real(dpc), intent(in) :: w
	end function Rc_M
	function Rc_a(w,a,dmn)
		use nrtype
		implicit none
		integer(i4b), intent(in) :: dmn
		complex(dpc), dimension(dmn) :: Rc_a
		complex(dpc), dimension(dmn), intent(in) :: a 
		real(dpc), intent(in) :: w
	end function Rc_a
end interface Rc

interface
	function Rr(w,dmn)
		use nrtype
		implicit none
		integer(i4b), intent(in) :: dmn
		real(dp), dimension(dmn,dmn) :: Rr
		real(dp), intent(in) :: w
	end function
end interface

interface
	subroutine SetANdiag_KS(a,Andiag)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in) :: a
		real(dp), dimension(:,:), intent(out) :: Andiag
	end subroutine
end interface

end module