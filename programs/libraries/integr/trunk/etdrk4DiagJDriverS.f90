subroutine etdrk4DiagJDriverS(ti,ai,Ji,h,tf,af,Jf,f0,f1,f2,f3,e,e2,Nplt,integrator,SetNlin,SetANdiag)

use nrtype
use nrutil, only: assert_eq
use ifc_integr !, only: etdrk4Diag

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
		complex(dpc), dimension(:,:), intent(out) :: Andiag
	end subroutine
end interface
interface
	subroutine integrator(a,J,h,aout,Jout,f0,f1,f2,f3,e,e2,SetNlin,SetAndiag)
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
! Simple driver for etdrk4Diag. The functions f0-1,e,e2 are precomputed
! in calling routine.


integer(i4b) :: d, dum, Nsteps, plrt, i,j 
complex(dpc), dimension(size(ai)) :: a
real(dp) :: t

d=2*(assert_eq(size(f3), size(a),'etdrk4Diag-a')-1)
dum=assert_eq(size(a), size(Ji,1), size(Jf,1), 'etdrk4Diag-a' )

t=ti
a=ai

Nsteps=nint((tf-ti)/h)	! Calculate number of steps
plrt=floor(Real(NSteps,dp)/Real(Nplt,dp)) ! Calculate after how many steps taken we should export values

if (allocated(tSt)) deallocate(tSt) !Clear out old stored variables if necessary.
if (allocated(aSt)) deallocate(aSt)
allocate(tSt(Nplt+1)) !Allocate storage for saved values.
allocate(aSt(Nplt+1,size(ai)))

j=1
tSt(j)=ti
aSt(1,:)=ai

Jf=Ji ! Change this.

do i=1,Nsteps
	call integrator(a,Jf,h,a,Jf,f0,f1,f2,f3,e,e2,SetNlin,SetANdiag)
	t=t+h
	if (mod(i,plrt) == 0) then ! export some value
		j=j+1
		tSt(j)=t
		aSt(j,:)=a  
	end if
end do

af=a ! return final point if required

end subroutine etdrk4DiagJDriverS