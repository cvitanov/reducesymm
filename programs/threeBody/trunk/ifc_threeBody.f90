module ifc_threeBody


interface
	subroutine CartToKSMcG_TBL(x,p_x,E,Q,P,r,Et,hyperR)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(in):: x,p_x
		real(dp), intent(in):: E 
		real(dp), dimension(:), intent(out) :: Q,P,r
		real(dp), intent(out)::Et, hyperR
	end subroutine
end interface

interface
	subroutine EOM_TBL(tau,f,v)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(in) :: f
		real(dp), intent(in) :: tau
		real(dp), dimension(:), intent(out) :: v
	end subroutine
end interface

interface
	subroutine KSMcGtoCart_TBL(Et,Q,P,hyperR,x,p_x,E)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(in) :: Q,P
		real(dp), intent(in)::  Et, hyperR
		real(dp), dimension(:), intent(out):: x,p_x
		real(dp), intent(out):: E 
	end subroutine
end interface


end module
