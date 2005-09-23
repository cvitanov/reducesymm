MODULE ifc_integr

use nrtype

! Contains all procedure declarations and allocation of temporary storage arrays.

! Storage of intermediate results to be communicated to the main programm or calling
! procedure.

! Needed by
real(dp), dimension(:), allocatable:: tSt
real(dp), dimension(:,:), allocatable :: ySt
complex(dpc), dimension(:,:), allocatable :: aSt



! Procedure declarations
interface
	subroutine derivsJ(x,y,J,dJds,MatVar)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: J
		REAL(DP), DIMENSION(:,:), INTENT(OUT) ::dJds
		interface
			function MatVar(x,y)
				USE nrtype
				IMPLICIT NONE
				REAL(DP), INTENT(IN) :: x
				REAL(DP), DIMENSION(:), INTENT(IN) :: y
				REAL(DP), DIMENSION(size(y),size(y)) :: MatVar
			end function MatVar
		end interface
	end subroutine
end interface

interface
	subroutine etdrk4Diag(a,h,aout,f0,f1,f2,f3,e,e2,SetNlin)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in):: a
		real(dp), intent(in) :: h
		real(dp), dimension(:),intent(in) :: f0,f1,f2,f3,e,e2
		complex(dpc), dimension(:), intent(out):: aout
		interface
			subroutine SetNlin(a,N_a)
			use nrtype
			implicit none
			complex(dpc), intent(in) :: a
			complex(dpc), intent(out) :: N_a
			end subroutine
		end interface
	end subroutine
end interface

interface
	subroutine etdrk4DiagDriverS(ti,ai,h,tf,af,f0,f1,f2,f3,e,e2,Nplt,SetNlin)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in):: ai
		real(dp), intent(in) :: ti,h,tf
		real(dp), dimension(:),intent(in) :: f0,f1,f2,f3,e,e2
		complex(dpc), dimension(:), intent(out) :: af
		integer(i4b), intent(in) :: Nplt
		interface
			subroutine SetNlin(a,N_a)
			use nrtype
			implicit none
			real(dpc), intent(in) :: a
			real(dpc), intent(out) :: N_a
			end subroutine
		end interface
	end subroutine
end interface

interface
	subroutine etdrk4DiagDriverS_a(ti,ai,h,tf,af,f0,f1,f2,f3,e,e2,Nplt,SetNlin)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in):: ai
		real(dp), intent(in) :: ti,h,tf
		real(dp), dimension(:),intent(in) :: f0,f1,f2,f3,e,e2
		complex(dpc), dimension(:), intent(out) :: af
		integer(i4b), intent(in) :: Nplt
		interface
			subroutine SetNlin(a,N_a)
			use nrtype
			implicit none
			real(dpc), intent(in) :: a
			real(dpc), intent(out) :: N_a
			end subroutine
		end interface
	end subroutine
end interface

interface
	subroutine etdrk4DiagJ(a,J,h,aout,Jout,f0,f1,f2,f3,e,e2,SetNlin,SetAndiag)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in):: a ! Initial point
		complex(dpc), dimension(:,:), intent(in):: J
		real(dp), intent(in) :: h ! Step size
		complex(dpc), dimension(:,:), intent(out):: Jout
		real(dp), dimension(:),intent(in) :: f0,f1,f2,f3,e,e2 ! Functions of the linear operator
		complex(dpc), dimension(:), intent(out):: aout ! Final point
		interface
			subroutine SetNlin(a,N_a)
			use nrtype
			implicit none
			complex(dpc), dimension(:), intent(in) :: a
			complex(dpc), dimension(:), intent(out) :: N_a
			end subroutine
		end interface
		interface
			subroutine SetAndiag(a,Andiag)
			use nrtype
			implicit none
			complex(dpc), dimension(:), intent(in) :: a
			complex(dpc), dimension(:,:), intent(out) :: Andiag
			end subroutine
		end interface
	end subroutine
end interface


interface
	subroutine etdrk4DiagJDriverS(ti,ai,Ji,h,tf,af,Jf,f0,f1,f2,f3,e,e2,Nplt,SetNlin,SetANdiag)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in):: ai ! Initial point
		complex(dpc), dimension(:,:), intent(in):: Ji 
		real(dp), intent(in) :: ti,h,tf ! initial time, stepsize, final time
		complex(dpc), dimension(:,:), intent(out):: Jf 
		real(dp), dimension(:),intent(in) :: f0,f1,f2,f3,e,e2 ! precomputed functions of the linear operator
		complex(dpc), dimension(:), intent(out) :: af ! Final point
		integer(i4b), intent(in) :: Nplt ! Number of intermediate points to be exported
		interface
			subroutine SetNlin(a,N_a)
			use nrtype
			implicit none
			complex(dpc), dimension(:), intent(in) :: a
			complex(dpc), dimension(:), intent(out) :: N_a
			end subroutine
		end interface
		interface
			subroutine SetANdiag(a,ANdiag)
			use nrtype
			implicit none
			complex(dpc), dimension(:), intent(in) :: a
			complex(dpc), dimension(:,:), intent(out) :: ANdiag
			end subroutine
		end interface
	end subroutine
end interface

interface
	subroutine etdrk4DiagPrefactors(Lin,h,R,M,f0,f1,f2,f3,e,e2)
		use nrtype
		implicit none
		real(dp), intent(in) :: Lin(:)
		real(dp), intent(in) :: h,R
		integer(i4b),intent(in) :: M
		real(dp), dimension(:),intent(out) :: f0,f1,f2,f3,e,e2
	end subroutine
end interface

interface
	subroutine eulerJ(x,y,J,Jout,h,MatVar)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), INTENT(IN) :: x,h
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: J
		REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Jout
		interface
			function MatVar(x,y)
				USE nrtype
				IMPLICIT NONE
				REAL(DP), INTENT(IN) :: x
				REAL(DP), DIMENSION(:), INTENT(IN) :: y
				REAL(DP), DIMENSION(size(y),size(y)) :: MatVar
			end function MatVar
		end interface
	end subroutine
end interface


interface
	subroutine integrPrk4(yi,Delta_x,qfP,yP,nsteps,nstepsP,nInters,sect,direction,derivs)
		USE nrtype 
		IMPLICIT NONE
		integer(i4b), intent(in) :: nsteps, nstepsP, nInters, sect
		real(dp), intent(in) ::  yi(:), Delta_x, qfP, direction
		real(dp), intent(out) :: yP(:,:)
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
				USE nrtype
				IMPLICIT NONE
				REAL(dp), INTENT(IN) :: x
				REAL(dp), DIMENSION(:), INTENT(IN) :: y
				REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx	
			END SUBROUTINE derivs
		END INTERFACE
	end subroutine
end interface

INTERFACE
	SUBROUTINE odeint(ystart,x1,x2,eps,h1,hmin,derivs,rkqs)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(INOUT) :: ystart
		REAL(DP), INTENT(IN) :: x1,x2,eps,h1,hmin
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
			REAL(DP), INTENT(IN) :: x
			REAL(DP), DIMENSION(:), INTENT(IN) :: y
			REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
	!BL
			SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
			USE nrtype
			REAL(DP), DIMENSION(:), INTENT(INOUT) :: y
			REAL(DP), DIMENSION(:), INTENT(IN) :: dydx,yscal
			REAL(DP), INTENT(INOUT) :: x
			REAL(DP), INTENT(IN) :: htry,eps
			REAL(DP), INTENT(OUT) :: hdid,hnext
				INTERFACE
				SUBROUTINE derivs(x,y,dydx)
					USE nrtype
					REAL(DP), INTENT(IN) :: x
					REAL(DP), DIMENSION(:), INTENT(IN) :: y
					REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
					END SUBROUTINE derivs
				END INTERFACE
			END SUBROUTINE rkqs
		END INTERFACE
	END SUBROUTINE odeint
END INTERFACE

interface
	Subroutine rk2J(x,y,dydx,h,yout,Ji,Jout,MatVar,derivs)
		USE nrtype 
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: y,dydx
		REAL(DP), INTENT(IN) :: x,h
		REAL(DP), DIMENSION(:), INTENT(OUT) :: yout
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: Ji
		REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Jout
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
	end subroutine
end interface

interface
	SUBROUTINE rk2Jdriver(xi,yi,xf,nsteps,y,Ji,Jout,MatVar,derivs)
		USE nrtype
		IMPLICIT none
		REAL(DP), INTENT(IN) :: xi,xf
		REAL(DP), DIMENSION(:), INTENT(IN) :: yi
		REAL(DP), DIMENSION(:), INTENT(OUT) :: y
		INTEGER(I4B), INTENT(IN) :: nsteps
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
	END SUBROUTINE
end interface

interface
	SUBROUTINE rk2Jsdriver(xi,yi,xf,nsteps,y,Ji,Jout,MatVar,derivs)
		USE nrtype
		IMPLICIT none
		REAL(DP), INTENT(IN) :: xi,xf
		REAL(DP), DIMENSION(:), INTENT(IN) :: yi
		REAL(DP), DIMENSION(:,:), INTENT(OUT) :: y
		INTEGER(I4B), INTENT(IN) :: nsteps
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
	END SUBROUTINE
end interface


INTERFACE
	SUBROUTINE rk4(x,y,dydx,h,yout,derivs)
		USE nrtype
		IMPLICIT NONE
		REAL(dp), DIMENSION(:), INTENT(IN) :: y,dydx
		REAL(dp), INTENT(IN) :: x,h
		REAL(dp), DIMENSION(:), INTENT(OUT) :: yout
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
				USE nrtype
				IMPLICIT NONE
				REAL(dp), INTENT(IN) :: x
				REAL(dp), DIMENSION(:), INTENT(IN) :: y
				REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs	
		END INTERFACE
	END SUBROUTINE
END INTERFACE


INTERFACE
	SUBROUTINE rk4driver(xi,yi,xf,nsteps,y,derivs)
		USE nrtype
		IMPLICIT none
		REAL(dp), INTENT(IN) :: xi,xf
		REAL(dp), DIMENSION(:), INTENT(IN) :: yi
		REAL(dp), DIMENSION(:,:), INTENT(OUT) :: y
		INTEGER(I4B), INTENT(IN) :: nsteps
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
				USE nrtype
				IMPLICIT NONE
				REAL(dp), INTENT(IN) :: x
				REAL(dp), DIMENSION(:), INTENT(IN) :: y
				REAL(dp), DIMENSION(:), INTENT(OUT) :: dydx	
			END SUBROUTINE derivs
		END INTERFACE
	END SUBROUTINE
END INTERFACE

interface
	Subroutine rk4J(x,y,dydx,h,yout,J,dJds,Jout,MatVar,derivs)
		USE nrtype 
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: y,dydx
		REAL(DP), INTENT(IN) :: x,h
		REAL(DP), DIMENSION(:), INTENT(OUT) :: yout
		REAL(DP), DIMENSION(:,:), INTENT(IN) :: J, dJds
		REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Jout
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
	end subroutine
end interface

interface
	SUBROUTINE rk4Jdriver(xi,yi,xf,nsteps,y,Ji,Jout,MatVar,derivs)
		USE nrtype
		IMPLICIT none
		REAL(DP), INTENT(IN) :: xi,xf
		REAL(DP), DIMENSION(:), INTENT(IN) :: yi
		REAL(DP), DIMENSION(:), INTENT(OUT) :: y
		INTEGER(I4B), INTENT(IN) :: nsteps
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
	END SUBROUTINE
end interface

INTERFACE
	SUBROUTINE rk4P(x,y,dydx,h,yout,derivs,p,sect)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: y,dydx
		REAL(DP), INTENT(IN) :: x,h
		REAL(DP), DIMENSION(:), INTENT(OUT) :: yout
		INTEGER(I4B), INTENT(IN) :: p,sect
		INTERFACE
			SUBROUTINE derivs(x,y,dydx,kappa)
				USE nrtype
				IMPLICIT NONE
				REAL(DP), INTENT(IN) :: x,kappa
				REAL(DP), DIMENSION(:), INTENT(IN) :: y
				REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx	
			END SUBROUTINE derivs
		END INTERFACE
	END SUBROUTINE rk4P
END INTERFACE

INTERFACE
	SUBROUTINE rk4Pdriver(xi,yi,xf,nsteps,y,derivs,p,sect)
		USE nrtype
		IMPLICIT none
		REAL(DP), INTENT(IN) :: xi,xf
		REAL(DP), DIMENSION(:), INTENT(IN) :: yi
		REAL(DP), DIMENSION(:,:), INTENT(OUT) :: y
		INTEGER(I4B), INTENT(IN) :: nsteps,p,sect
		INTERFACE
			SUBROUTINE derivs(x,y,dydx,kappa)
				USE nrtype
				IMPLICIT NONE
				REAL(DP), INTENT(IN) :: x,kappa
				REAL(DP), DIMENSION(:), INTENT(IN) :: y
				REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx	
			END SUBROUTINE derivs
		END INTERFACE
	END SUBROUTINE rk4Pdriver
END INTERFACE

INTERFACE
	SUBROUTINE rkck(y,dydx,x,h,yout,yerr,derivs)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(IN) :: y,dydx
		REAL(DP), INTENT(IN) :: x,h
		REAL(DP), DIMENSION(:), INTENT(OUT) :: yout,yerr
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
			REAL(DP), INTENT(IN) :: x
			REAL(DP), DIMENSION(:), INTENT(IN) :: y
			REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
	END SUBROUTINE rkck
END INTERFACE

INTERFACE
	SUBROUTINE rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
		USE nrtype
		REAL(DP), DIMENSION(:), INTENT(INOUT) :: y
		REAL(DP), DIMENSION(:), INTENT(IN) :: dydx,yscal
		REAL(DP), INTENT(INOUT) :: x
		REAL(DP), INTENT(IN) :: htry,eps
		REAL(DP), INTENT(OUT) :: hdid,hnext
		INTERFACE
			SUBROUTINE derivs(x,y,dydx)
			USE nrtype
			REAL(DP), INTENT(IN) :: x
			REAL(DP), DIMENSION(:), INTENT(IN) :: y
			REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx
			END SUBROUTINE derivs
		END INTERFACE
	END SUBROUTINE rkqs
END INTERFACE


END MODULE