module ifc_newton

interface
	subroutine NewtonPO(yi,T,a,tol,maxIter,nsteps,derivs,MatVar,conv,J)
		use nrtype
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
	end subroutine
end interface

interface
	subroutine setLHM(J,v,a,LHM)
		use nrtype
		implicit none
		real(dp), dimension(:,:), intent(in) :: J
		real(dp), dimension(:), intent(in) :: a, v
		real(dp), dimension(:,:), intent(out) :: LHM 
	end subroutine
end interface

interface
	subroutine setRHS(diff,RHS)
		use nrtype
		implicit none
		real(dp),  dimension(:), intent(in) :: diff
		real(dp),  dimension(:), intent(out) :: RHS
	end subroutine
end interface

end module