	SUBROUTINE TBrkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	USE ifc_integr, ONLY : rkck
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: y
	REAL(DP), DIMENSION(:), INTENT(IN) :: dydx,yscal
	REAL(DP), INTENT(INOUT) :: x
	REAL(DP), INTENT(IN) :: htry,eps
	REAL(DP), INTENT(OUT) :: hdid,hnext
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B) :: ndum
	REAL(DP) :: errmax,h,htemp,xnew
	REAL(DP), DIMENSION(size(y)) :: yerr,ytemp
	REAL(DP), PARAMETER :: SAFETY=0.9_dp,PGROW=-0.2_dp,PSHRNK=-0.25_dp,&
		ERRCON=1.89e-4
	ndum=assert_eq(size(y),size(dydx),size(yscal),'rkqs')
	h=htry
	do
		call rkck(y,dydx,x,h,ytemp,yerr,derivs)
		if ( (Abs(1/(y(1)**2+y(2)**2)) - Abs(yerr(5)) ) <= 0 ) call nrerror('Energy trancution error larger than e-e interaction term.')
		errmax=maxval(abs(yerr(:)/yscal(:)))/eps
		if (errmax <= 1.0) exit
		htemp=SAFETY*h*(errmax**PSHRNK)
		h=sign(max(abs(htemp),0.1_dp*abs(h)),h)
		xnew=x+h
		if (xnew == x) call nrerror('stepsize underflow in TBrkqs')
	end do
	if (errmax > ERRCON) then
		hnext=SAFETY*h*(errmax**PGROW)
	else
		hnext=5.0_dp*h
	end if
	hdid=h
	x=x+h
	y(:)=ytemp(:)
	END SUBROUTINE TBrkqs
