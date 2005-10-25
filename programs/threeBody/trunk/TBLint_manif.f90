MODULE TBL_path_manif
	USE nrtype
	INTEGER(I4B) :: nok,nbad,kount
	LOGICAL(LGT), SAVE :: save_steps=.false.
	REAL(DP) :: dxsav
	REAL(DP), DIMENSION(:), POINTER :: xp
	REAL(DP), DIMENSION(:,:), POINTER :: yp
END MODULE TBL_path_manif

	SUBROUTINE TBLint_manif(ystart,x1,x2,eps,h1,hmin,derivs,rkqs,pChange)
	USE nrtype; USE nrutil, ONLY : nrerror,reallocate
	USE TBL_path_manif
	USE ifc_threeBody
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: ystart
	REAL(DP), INTENT(IN) :: x1,x2,eps,h1,hmin
	logical(lgt), intent(out) :: pChange
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
!BL
		SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
		USE nrtype
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
		END SUBROUTINE rkqs
	END INTERFACE
	REAL(DP), PARAMETER :: TINY=1.0e-30_dp, SAFEF=1.05_dp, xasym = 200.0_dp
	INTEGER(I4B), PARAMETER :: MAXSTP=100000000
	INTEGER(I4B) :: nstp
	REAL(DP) :: h,hdid,hnext,x,xsav
	REAL(DP), DIMENSION(size(ystart)) :: dydx,y,yscal
	real(dp), dimension(2) :: xx,p_x,xxstart
	real(dp) :: E

	pChange=.false.
	 
	call KSMcGtoCart_TBL(ystart(5),ystart(1:2),ystart(3:4),ystart(6),xxstart,p_x,E)
	x=x1
	h=sign(h1,x2-x1)
	nok=0
	nbad=0
	kount=0
	y(:)=ystart(:)
	if (save_steps) then
		xsav=x-2.0_dp*dxsav
		nullify(xp,yp)
		allocate(xp(256))
		allocate(yp(size(ystart),size(xp)))
	end if
	do nstp=1,MAXSTP
		call derivs(x,y,dydx)
		yscal(:)=abs((h)*dydx(:))+hmin
		if (save_steps .and. (abs(x-xsav) > abs(dxsav)))  then 
			call save_a_step
		end if
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x
		call rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs)
		if (hdid == h) then
			nok=nok+1
		else
			nbad=nbad+1
		end if
		call KSMcGtoCart_TBL(y(5),y(1:2),y(3:4),y(6),xx,p_x,E)
!		print *,x,xx
		if ( p_x(1) < 0.0_dp ) then
			pChange=.true.
			ystart(:)=y(:)
			if (save_steps) call save_a_step
			RETURN	
		end if
		if (  ( (x-x2)*(x2-x1) >= 0.0 ) .or. ( MaxVal(Abs(xx)) >=  xasym ) ) then ! then  !! Edited for TBL
!			print *,x
			ystart(:)=y(:)
			if (save_steps) call save_a_step
			RETURN
		end if
		if (abs(hnext) < hmin) then 
			call nrerror('stepsize smaller than minimum in odeint')
		end if
		h=hnext
	end do
	call nrerror('too many steps in odeint')
	CONTAINS
!BL
	SUBROUTINE save_a_step
	kount=kount+1
	if (kount > size(xp)) then
		xp=>reallocate(xp,2*size(xp))
		yp=>reallocate(yp,size(yp,1),size(xp))
	end if
	xp(kount)=x
	yp(:,kount)=y(:)
	xsav=x
	END SUBROUTINE save_a_step
	END SUBROUTINE TBLint_manif
