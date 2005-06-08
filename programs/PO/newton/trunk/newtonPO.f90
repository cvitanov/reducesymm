subroutine NewtonPO(y,T,a,tol,maxIter,nsteps,derivs,MatVar,conv)

real(dp) :: intent(inout):: yi(:)
real(dp) :: intent(inout) :: T
real(dp) :: intent(in) :: a(:)
integer(i4b) :: intent(in) :: sect, maxIter,nsteps
integer(i4b) :: intent(out) :: conv
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
real(dp) :: yi(size(yi)), xi, xf, v
integer(i4b) :: i, ndum, ic
real(dp) :: a(size(y)-1),  indx(size(y))
real(dp) :: LHM(size(y),size(y)), RHS(size(y))

ndum=assert_eq(size(y)-1,size(a),'NewtonPO')

xi=0.0_dp

conv=0

do i=1,maxIter
	call rk4Jdriver(xi,yi,T,nsteps,y,Ji,Jout,MatVar,derivs)
	if (max(abs(y-yi))) then
		conv=1
		exit
	end if 
	call derivs(y(size(y)),y,v)
	call setLHM(J,v(1:size(v)-1),a,LHM)
	call setRHS(yi(1:size(yi)-1),y(1:size(y)-1),RHS)
	call ludcmp(LHM,indx,ic)
	call lubksb(LHM,indx,RHS)
	yi(1:size(yi)-1) = y(1:size(y)-1) + RHS(y(1:size(RHS)-1) ! update y and T
	yi(size(yi)) = 0.0_dp
	T = T - RHS(size(RHS))
end do

end subroutine