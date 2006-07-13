subroutine newtonPOetdrk4_g(ai,T,q,tol,maxIter,h,f0,f1,f2,f3,e,e2,SetLin,SetNlin,SetANdiag,conv,J)

use nrtype
use nrutil, only: assert_eq
use ifc_util
use ifc_integr
use ifc_po, only: setLHM, setRHS
use la_precision, only: wp => dp
use f95_lapack, only: LA_GESV, LA_GEEV

implicit none

complex(dpc), intent(inout):: ai(:)
real(dp), intent(inout) :: T		! Guess period in, computed period out
complex(dpc), intent(in) :: q(:) 	! Vector normal to Poincare hyperplane
real(dp), intent(in) :: tol		
integer(i4b), intent(in) :: maxIter	! Maximum number of iterations to be attempted
real(dp), intent(in) :: h
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
complex(dpc), dimension(size(ai)) :: af,diff, v
complex(dpc), dimension(2*(size(ai)+1)) :: RHS 
real(dp) :: ddum, ti, mx 
integer(i4b) :: i, ndum, Nplt=1, k,m
complex(dpc), dimension(2*(size(ai)+1),2*(size(ai)+1)) :: LHM, LHMdum
real(dp), dimension(size(ai)) :: Lin
complex(dpc), dimension(size(ai)) :: N_a
complex(dpc), dimension(size(J,1),size(J,2)):: Jdum
complex(dpc), dimension(size(J,1)):: Jeig
complex(dpc), dimension(size(J,1)+1):: LHMeig

ndum=assert_eq(size(ai),size(q)/2-1,size(J,1),size(J,2),'NewtonPO')

ti=0.0_dp

conv=0

do i=1,maxIter
	print *,i
	J = UnitMatrix(size(J,1))
	call etdrk4DiagJDriverS(ti,ai,J,h,T,af,J,f0,f1,f2,f3,e,e2,Nplt,etdrk4DiagJ_g,SetNlin,SetANdiag)
 	Jdum=J
	call la_geev(Jdum,Jeig)
 	print *,"*********************Jeig"
 	print *,Jeig
 	print *,"*********************"
	diff = af-ai
	mx=maxval(Abs(diff))
	print *,"Max",mx
	print *,"a2i",ai(3),"a2f",af(3),"diff",af(3)-ai(3)
	if (mx < tol) then ! Perharps need an extra confition 
		conv=1
		exit
	end if
	call setLin(Lin)
	call setNlin(af,N_a)
	v = Lin*af+N_a
	call setLHM(af,J,v,q,LHM)
	call setRHS(diff,RHS)
	call la_gesv(LHM,RHS)
	ai = ai + RHS(1:size(ai)) ! update y, T
	T = T - RHS(size(ai)+1)
end do

end subroutine newtonPOetdrk4_g



subroutine newtonPOetdrk4_c0(ai,T,q,tol,maxIter,h,f0,f1,f2,f3,e,e2,SetLin,SetNlin,SetANdiag,conv,J,c0)

use nrtype
use nrutil, only: assert_eq
use ifc_util
use ifc_integr
use ifc_po, only: setLHM, setRHS
use la_precision, only: wp => dp
use f95_lapack, only: LA_GESV, LA_GEEV

implicit none

complex(dpc), intent(inout):: ai(:)
real(dp), intent(inout) :: T		! Guess period in, computed period out
complex(dpc), intent(in) :: q(:) 	! Vector normal to Poincare hyperplane
real(dp), intent(in) :: tol		
integer(i4b), intent(in) :: maxIter	! Maximum number of iterations to be attempted
real(dp), intent(in) :: h
real(dp), dimension(:),intent(in) :: f0,f1,f2,f3,e,e2 ! precomputed functions of the linear operator
integer(i4b), intent(out) :: conv	! If converged conv -> 1, if not conv -> 0
complex(dpc), intent(out) :: J(:,:)
integer(i4b), intent(in), optional:: c0 !
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
complex(dpc), dimension(size(ai)) :: af,diff, v
complex(dpc), dimension(2*(size(J,1)+1)) :: RHS 
real(dp) :: ddum, ti, mx 
integer(i4b) :: i, ndum, Nplt=1, k,m
complex(dpc), dimension(2*(size(J,1)+1),2*(size(J,2)+1)) :: LHM, LHMdum
real(dp), dimension(size(ai)) :: Lin
complex(dpc), dimension(size(ai)) :: N_a
complex(dpc), dimension(size(J,1),size(J,2)):: Jdum
complex(dpc), dimension(size(J,1)):: Jeig
complex(dpc), dimension(size(J,1)+1):: LHMeig

ndum=assert_eq(size(ai),size(q)/2,size(J,1)+1,size(J,2)+1,'NewtonPO')

ti=0.0_dp

conv=0

do i=1,maxIter
	print *,i
	J = UnitMatrix(size(J,1))
	call etdrk4DiagJDriverS(ti,ai,J,h,T,af,J,f0,f1,f2,f3,e,e2,Nplt,etdrk4DiagJ_c0,SetNlin,SetANdiag,c0=1)
print *,"after"
 	Jdum=J
	call la_geev(Jdum,Jeig)
 	print *,"*********************Jeig"
 	print *,Jeig
 	print *,"*********************"
	diff = af-ai
	mx=maxval(Abs(diff))
	print *,"Max",mx
	print *,"a2i",ai(3),"a2f",af(3),"diff",af(3)-ai(3)
	if (mx < tol) then ! Perharps need an extra confition 
		conv=1
		exit
	end if
	call setLin(Lin)
	call setNlin(af,N_a)
	v = Lin*af+N_a
	call setLHM(af(2:size(af)),J,v(2:size(af)),q,LHM)
	call setRHS(diff(2:size(af)),RHS)
	call la_gesv(LHM,RHS)
	print *,"dP",RHS(2)
	ai(2:size(ai))= ai(2:size(ai)) + RHS(1:size(J,1)) ! update y, T
	print *,"deltaT",RHS(size(RHS))
	print *,"deltaa2",RHS(2)
	T = T - RHS(size(RHS))
	print *,"T",T
end do

end subroutine newtonPOetdrk4_c0