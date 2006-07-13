module ifc_threeBody


interface
	subroutine CartToKSMcG_TBL(x,p_x,E,Q,P,r,Et,hyperR)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(in):: x,p_x
		real(dp), intent(in):: E 
		real(dp), dimension(:), intent(out) :: Q,P,r
		real(dp), intent(out)::Et, hyperR
	end subroutine
end interface

interface
	subroutine findIC(xi,p_xi,E,Z)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(inout) :: xi
		real(dp), dimension(:), intent(in) :: p_xi
		real(dp), intent(in) :: E
		real(dp), intent(in) :: Z
	end subroutine
end interface

interface findSep
	subroutine findSep_TBL(E,Z,xo,xmax,pUp,pDown,tol,Npoints,MaxAttempts,eps,h1,hmin,ti,taui,tauf,EOM_TBL)
		use nrtype
		implicit none
		real(dp), intent(in) :: E, Z, xo,xmax
		real(dp), dimension(:), intent(inout) ::  pUp, pDown
		real(dp), intent(in) :: tol 
		integer(i4b), intent(in) :: Npoints, MaxAttempts
		real(dp), intent(in) :: eps, h1, hmin
		real(dp), intent(in) :: ti, taui, tauf
		interface
			subroutine EOM_TBL(tau,f,v)
				use nrtype
				implicit none
				real(dp), dimension(:), intent(in) :: f
				real(dp), intent(in) :: tau
				real(dp), dimension(:), intent(out) :: v
			end subroutine
		end interface
	end subroutine
end interface

interface
	subroutine KSMcGtoCart_TBL(Et,Q,P,hyperR,x,p_x,E)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(in) :: Q,P
		real(dp), intent(in)::  Et, hyperR
		real(dp), dimension(:), intent(out):: x,p_x
		real(dp), intent(out):: E 
	end subroutine
end interface

interface 
	subroutine TBLint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(INOUT) :: ystart
		REAL(DP), INTENT(IN) :: x1,x2,eps,h1,hmin
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
	end subroutine
end interface

interface 
	subroutine TBLint_sep(ystart,x1,x2,eps,h1,hmin,derivs,rkqs,findMan)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(INOUT) :: ystart
		REAL(DP), INTENT(IN) :: x1,x2,eps,h1,hmin
		logical(lgt), intent(out) :: findMan
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
	end subroutine
end interface


interface
	SUBROUTINE TBrkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
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
	END SUBROUTINE TBrkqs
end interface

interface
	subroutine vectorize(Q,P,Et,hyperR,ti,fi)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(in):: Q,P
		real(dp), intent(in):: Et,hyperR,ti
		real(dp), dimension(:), intent(out):: fi 
	end subroutine vectorize
end interface

end module
