subroutine newtonPOetdrk4(ai,T,kappa,q,tol,maxIter,h,diagk,f0,f1,f2,f3,e,e2,SetLin,SetNlin,SetANdiag,conv,J)

use nrtype
use nrutil, only: assert_eq
use ifc_util, only: UnitMatrix
use ifc_lu, only: lubksb, ludcmp
use ifc_integr, only: etdrk4DiagJDriverS, etdrk4DiagJ
use ifc_rpo, only: setLHM, setRHS
use la_precision, only: wp => dp
use f95_lapack, only: la_gesv

implicit none

complex(dpc), intent(inout):: ai(:)
real(dp), intent(inout) :: T		! Guess period in, computed period out
real(dp), intent(inout) :: kappa	! Guess kappa in, computed kappa out
complex(dpc), intent(in) :: q(:) 	! Vector normal to Poincare hyperplane
real(dp), intent(in) :: tol		
integer(i4b), intent(in) :: maxIter	! Maximum number of iterations to be attempted
real(dp), intent(in) :: h
real(dp), dimension(:), intent(in) :: diagk	! Diag[ik/tilde{L}]	
real(dp), dimension(:),intent(in) :: f0,f1,f2,f3,e,e2 ! precomputed functions of the linear operator
integer(i4b), intent(out) :: conv	! If converged conv -> 1, if not conv -> 0
complex(dpc), intent(out) :: J(:,:)
interface
        subroutine SetLin(Lin)
                use nrtype
                implicit none
                real(dp), dimension(:), intent(out) :: Lin
        end subroutine
end interface
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
!
!
complex(dpc), dimension(size(ai)) :: af,diff, v, RHS 
real(dp) :: ddum, ti, tf, mx 
integer(i4b) :: i, ndum, Nplt=1
complex(dpc), dimension(size(ai),size(ai)) :: LHM
real(dp), dimension(size(ai)) :: Lin
complex(dpc), dimension(size(ai)) :: N_a, R

ndum=assert_eq(size(ai),size(q),size(J,1),size(J,2),'NewtonPO')

ti=0.0_dp

conv=0

do i=1,maxIter
	J = UnitMatrix(size(J,1))
	call etdrk4DiagJDriverS(ti,ai,J,h,tf,af,J,f0,f1,f2,f3,e,e2,Nplt,etdrk4DiagJ,SetNlin,SetANdiag)
	diff = R*af(1:size(af)-1)-ai(1:size(ai)-1)
	mx=maxval(Abs(diff))
	print *,mx
	if (mx < tol) then
		conv=1
		exit
	end if
	call setLin(Lin)
	call setNlin(af,N_a)
	v = Lin*af+N_a
	R = exp(kappa*diagk) 
	call setLHM(af,J,v,q,R,diagk,LHM)
	call setRHS(diff,RHS)
	call la_gesv(LHM,RHS)
	ai = ai + RHS(size(ai)) ! update y, T and kappa
	T = T + RHS(size(ai)+1)
	kappa = kappa + RHS(size(ai)+2) 
end do

end subroutine