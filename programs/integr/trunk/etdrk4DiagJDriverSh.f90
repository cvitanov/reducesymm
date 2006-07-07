subroutine etdrk4DiagJDriverSh(ti,ai,Ji,h,Nsteps,tf,af,Jf,f0,f1,f2,f3,e,e2,f0dum,f1dum,f2dum,f3dum,edum,e2dum,Nplt,integrator,integratorJ,SetNlin,SetANdiag)

use nrtype
use nrutil, only: assert_eq
use ifc_integr

implicit none

complex(dpc), dimension(:), intent(in):: ai ! Initial point
real(dp), dimension(:,:), intent(in):: Ji 
real(dp), intent(in) :: ti,h,tf ! initial time, stepsize, final time
real(dp), dimension(:,:), intent(out):: Jf 
real(dp), dimension(:),intent(in) :: f0,f1,f2,f3,e,e2 ! precomputed functions of the linear operator
real(dp), dimension(:),intent(in) :: f0dum,f1dum,f2dum,f3dum,edum,e2dum ! precomputed functions of the linear operator
complex(dpc), dimension(:), intent(out) :: af ! Final point
integer(i4b), intent(in) :: Nplt,Nsteps ! Number of intermediate points to be exported
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
		real(dp), dimension(:,:), intent(out) :: Andiag
	end subroutine
end interface
interface
	subroutine integrator(a,h,aout,f0,f1,f2,f3,e,e2,SetNlin)
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
			complex(dpc), dimension(:), intent(in) :: a
			complex(dpc), dimension(:), intent(out) :: N_a
			end subroutine
		end interface
	end subroutine
end interface
interface
	subroutine integratorJ(a,adum,af,J,h,Jout,f0,f1,f2,f3,e,e2,SetNlin,SetAndiag)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in):: a,adum,af ! Initial, middle and final point
		real(dp), dimension(:,:), intent(in):: J
		real(dp), intent(in) :: h ! Step size
		real(dp), dimension(:,:), intent(out):: Jout
		real(dp), dimension(:),intent(in) :: f0,f1,f2,f3,e,e2 ! Functions of the linear operator
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
			real(dp), dimension(:,:), intent(out) :: Andiag
			end subroutine
		end interface
	end subroutine
end interface
! Driver to compute the Jacobian using etdrk4Diag. The functions f0-1*,e*,e2* (dum for half the stepsize) 
! are precomputed in calling routine.


integer(i4b) :: d, dum, plrt, i,j 
complex(dpc), dimension(size(ai)) :: a,adum,ao
real(dp) :: t, h2

d=2*(assert_eq(size(f3), size(a),'etdrk4DiagJDriver-a')-1)
dum=assert_eq(d, size(Ji,1), size(Jf,1), 'etdrk4DiagJDriver-J' )

t=ti
a=ai

plrt=ceiling(Real(NSteps,dp)/Real(Nplt,dp)) ! Calculate after how many steps taken we should export values

if (allocated(tSt)) deallocate(tSt) !Clear out old stored variables if necessary.
if (allocated(aSt)) deallocate(aSt)
allocate(tSt(Nplt+1)) !Allocate storage for saved values.
allocate(aSt(Nplt+1,size(ai)))

tSt=0.0_dp
aSt=(0.0,0.0)

!print *,"alloc",allocated(tSt),allocated(aSt)

j=1
tSt(j)=ti
aSt(1,:)=ai

Jf=Ji ! Change this.

h2=h/2.0_dp

do i=1,Nsteps
	call integrator(a,h2,adum,f0dum,f1dum,f2dum,f3dum,edum,e2dum,SetNlin)
	call integrator(a,h,ao,f0,f1,f2,f3,e,e2,SetNlin)
	call integratorJ(a,adum,ao,Jf,h,Jf,f0,f1,f2,f3,e,e2,SetNlin,SetANdiag)
	t=ti+i*h
	a=ao
	if (mod(i,plrt) == 0) then ! export some value
		j=j+1
		tSt(j)=t
		aSt(j,:)=a
	end if
end do

af=a ! return final point if required

end subroutine etdrk4DiagJDriverSh
