subroutine integrP(yi,Delta_x,qfP,yP,nsteps,nstepsP,nInters,sect,derivs)

USE nrtype; USE ifc_integr

IMPLICIT NONE

integer(i4b), intent(in) :: nsteps, nstepsP, nInters, sect
real(dp), intent(in) ::  yi(:), Delta_x, qfP
real(dp), intent(out) :: yP(:,:)

REAL(DP), DIMENSION(:,:), ALLOCATABLE :: y, y_poinc
real(dp) :: xi
INTEGER(I4B) :: d, p, iI = 0 
INTEGER(I4B) :: idum, ndum, i
INTERFACE
	SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), INTENT(IN) :: x
		REAL(dp), DIMENSION(:), INTENT(IN) :: y
		REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE derivs
END INTERFACE



ndum = assert_eq(size(yP,1), nstepsP, "integrP: nstepsP")
d = assert_eq(size(yi),size(yP,2, "integrP: d-1")) - 1 

allocate(y(nsteps,d+1))
allocate(y_poinc(nstepsP,d+1))

p=0

y=0.0_dp

do i=1,d+1
	y(1,i)=yi(i)
end do

xi=y(1,d+1)

do while( (iI<nInters) )
	call rk4driver(y(1,d+1),y(1,:),y(1,d+1)+Delta_x,nsteps,y,derivs,p,sect)
	if  ( ( y(nsteps,sect) > qfP ) .AND. ( y(1,sect) < qfP  ) ) THEN ! Intersected Poincare, refine.
		p=1
		call rk4driver(y(1,sect),y(1,:),qfP,nstepsP,y_poinc,derivs,p,sect)
		iI=iI+1
		p=0
		yP(iI) = y_poinc(nstepsP,:)
		y=0.0_dp
		y(1,:)=y_poinc(nstepsP,:)
		y(1,sect)=qfP
		y_poinc=0.0_dp
	else 		! Forget history and continue.
		y(1,:)=y(nsteps,:)
		y(2:nsteps,:)=0.0_dp 
	end if
end do



end subroutine