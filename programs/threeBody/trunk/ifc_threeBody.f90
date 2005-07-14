module ifc_threeBody


interface
	subroutine CartToKSMcG_TBL(x,p_x,dt,E,Q,P,r,dtau,Et,hyperR)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(in):: x,p_x
		real(dp), intent(in):: dt, E 
		real(dp), dimension(:), intent(out) :: Q,P,r
		real(dp), intent(out):: dtau, Et, hyperR
	end subroutine
end interface

interface
	subroutine EOM_TBL(Q,P,r,E,hyperR,v)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(in) :: Q,P,r
		real(dp), intent(in) :: E,hyperR
		real(dp), dimension(:), intent(out) :: v
	end subroutine
end interface

interface
	subroutine KSMcGtoCart_TBL(Et,Q,P,hyperR,x,p_x,E)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(in) :: Q,P,r
		real(dp), intent(in)::  Et, hyperR
		real(dp), dimension(:), intent(out):: x,p_x
		real(dp), intent(out):: E 
	end subroutine
end interface


end module