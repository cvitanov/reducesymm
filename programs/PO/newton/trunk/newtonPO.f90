subroutine NewtonPO(yi,T,a,tol,maxIter,nsteps,derivs,MatVar,conv)

use nrtype; use nrutil, only: assert_eq
use ifc_util, only: UnitMatrix

implicit none

real(dp), intent(inout):: yi(:)
real(dp), intent(inout) :: T
real(dp), intent(in) :: a(:), tol
integer(i4b), intent(in) :: maxIter,nsteps
integer(i4b), intent(out) :: conv
interface
	subroutine derivs(x,y,dydx)
		use nrtype
		implicit none
		real(dp), intent(in) :: x
		real(dp), dimension(:), intent(in) :: y
		real(dp), dimension(:), intent(out) :: dydx	
	end subroutine derivs
end interface
interface
	function MatVar(x,y)
		use nrtype
		implicit none
		real(dp), intent(in) :: x
		real(dp), dimension(:), intent(in) :: y
		real(dp), dimension(size(y),size(y)) :: MatVar
	end function MatVar
end interface
!
!
real(dp) :: y(size(yi)),diff(size(yi)), xi, xf, mx, mn, v(size(yi))
real(dp), dimension(size(yi)-1,size(yi)-1) :: J
integer(i4b) :: i, ndum, ic
real(dp) :: indx(size(yi))
real(dp) :: LHM(size(yi),size(yi)), RHS(size(yi))

ndum=assert_eq(size(y)-1,size(a),'NewtonPO')

xi=0.0_dp

conv=0

do i=1,maxIter
	J = UnitMatrix(size(J,1))
	call rk4Jdriver(xi,yi,T,nsteps,y,J,J,MatVar,derivs)
	diff = y-yi
	mx=maxval(diff)
	mn=minval(diff)
	mx=max(abs(mx),abs(mn))
	if (mx < tol) then
		conv=1
		exit
	end if 
	call derivs(y(size(y)),y,v)
	call setLHM(J,v(1:size(v)-1),a,LHM)
	call setRHS(diff(1:size(diff-1)),RHS) 
	call ludcmp(LHM,indx,ic)
	call lubksb(LHM,indx,RHS)
	yi(1:size(yi)-1) = y(1:size(y)-1) + RHS(1:size(RHS)-1) ! update y and T
	yi(size(yi)) = 0.0_dp
	T = T - RHS(size(RHS))
end do

end subroutine