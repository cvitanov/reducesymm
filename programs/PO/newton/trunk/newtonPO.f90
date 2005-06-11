subroutine NewtonPO(yi,T,a,tol,maxIter,nsteps,derivs,MatVar,conv,J)

use nrtype
use nrutil, only: assert_eq
use ifc_util, only: UnitMatrix
use ifc_lu, only: lubksb, ludcmp
use ifc_integr, only: rk4Jdriver
use ifc_newton, only: setLHM, setRHS

implicit none

real(dp), intent(inout):: yi(:)
real(dp), intent(inout) :: T
real(dp), intent(in) :: a(:), tol
integer(i4b), intent(in) :: maxIter,nsteps
integer(i4b), intent(out) :: conv
real(dp), intent(out) :: J(:,:)
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
real(dp) :: y(size(yi)),diff(size(yi)-1), ddum, xi, xf, mx, mn, v(size(yi))
integer(i4b) :: i, ndum,  indx(size(yi))
real(dp) :: LHM(size(yi),size(yi)), RHS(size(yi))

ndum=assert_eq(size(y)-1,size(a),size(J,1),size(J,2),'NewtonPO')

xi=0.0_dp

conv=0

do i=1,maxIter
	J = UnitMatrix(size(J,1))
	call rk4Jdriver(xi,yi,T,nsteps,y,J,J,MatVar,derivs)
	diff = y(1:size(y)-1)-yi(1:size(yi)-1)
	mx=maxval(diff)
	mn=minval(diff)
	mx=max(abs(mx),abs(mn))
	print *,mx
	if (mx < tol) then
		conv=1
		exit
	end if
	call derivs(y(size(y)),y,v)
	call setLHM(J,v(1:size(v)-1),a,LHM)
	call setRHS(diff,RHS)
	call ludcmp(LHM,indx,ddum)
	call lubksb(LHM,indx,RHS)
	yi(1:size(yi)-1) = yi(1:size(y)-1) + RHS(1:size(RHS)-1) ! update y and T
	yi(size(yi)) = 0.0_dp
	T = T - RHS(size(RHS))
end do

if (conv == 1)	then
	yi(size(yi)) = 0.0_dp
	J = UnitMatrix(size(J,1))
	call rk4Jdriver(xi,yi,T,nsteps,y,J,J,MatVar,derivs)
end if

end subroutine