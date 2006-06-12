SUBROUTINE rk4JdriverS(xi,yi,xf,Nsteps,Nplt,y,Ji,Jout,MatVar,derivs)

USE nrtype ; USE ifc_integr
USE nrutil, ONLY : assert_eq

implicit none

REAL(DP), INTENT(IN) :: xi,xf
REAL(DP), DIMENSION(:), INTENT(IN) :: yi
REAL(DP), DIMENSION(:), INTENT(OUT) :: y
INTEGER(I4B), INTENT(IN) :: Nsteps,Nplt
REAL(DP), DIMENSION(:,:), INTENT(IN) :: Ji
real(dp), dimension(:,:), intent(out) :: Jout
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
!
!

REAL(DP) :: h, v(size(yi)), dJds(size(Ji,1),size(Ji,2))
INTEGER(I4B) :: i,j=1, mdum,plrt

print *,"yi",yi

mdum=assert_eq(size(y),size(yi),'rk2Jdriver: y dim2')

h=abs(xf-xi)/Nsteps

plrt=ceiling(Real(NSteps,dp)/Real(Nplt,dp)) ! Calculate after how many steps taken we should export values
print *,"t-rk4J",Nsteps,plrt,Nsteps*h,xf,h

y(:)=yi

if (allocated(tSt)) deallocate(tSt) !Clear out old stored variables if necessary.
if (allocated(ySt)) deallocate(ySt)
allocate(tSt(Nplt)) !Allocate storage for saved values.
allocate(ySt(Nplt,size(yi)-1))
tSt=0.0_dp
ySt=0.0_dp

Jout=Ji

DO i=2,nsteps+1
	call derivs(y(size(y)),y,v)
	call derivsJ(y(size(y)),y(1:size(y)-1), Jout, dJds, MatVar)
	call rk4J(y(size(y)),y,v,h,y,Jout,dJds,Jout,MatVar,derivs)
	if (mod(i,plrt) == 0) then ! export some value
!		print *,i,j
		tSt(j)=y(size(y))
		ySt(j,:)=y(1:size(y)-1)
		j=j+1
	end if
END DO



END SUBROUTINE