	FUNCTION rtbisB(func,x1,x2,xacc)
	USE nrtype; USE nrutil, ONLY : nrerror
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x1,x2,xacc
	REAL(DP) :: rtbis
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		logical(lgt) :: func
		END FUNCTION func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXIT=40
	INTEGER(I4B) :: j
	REAL(DP) :: dx,xmid
	logical(lgt) :: f,fmid, negative, positive
	fmid=func(x2)
	f=func(x1)
	if ( f .xor. fmid ) call nrerror('rtbis: root must be bracketed')
	negative=f
	positive=fmid
	rtbis=x1
	dx=x2-x1
	do j=1,MAXIT
		dx=dx*0.5_DP
		xmid=rtbis+dx
		fmid=func(xmid)
		if ( .not. (fmid .xor. 0.0) ) rtbis=xmid
		if (abs(dx) < xacc) RETURN
	end do
	call nrerror('rtbis: too many bisections')
	END FUNCTION rtbisB
