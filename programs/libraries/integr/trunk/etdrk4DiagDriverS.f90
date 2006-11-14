subroutine etdrk4DiagDriverS(ti,ai,Nsteps,tf,af,f0,f1,f2,f3,e,e2,Nplt,SetNlin)

use nrtype
use nrutil, only: assert_eq
use ifc_integr

implicit none

complex(dpc), dimension(:), intent(in):: ai ! Initial point
real(dp), intent(in) :: ti,tf ! initial time, stepsize, final time
real(dp), dimension(:),intent(in) :: f0,f1,f2,f3,e,e2 ! precomputed functions of the linear operator
complex(dpc), dimension(:), intent(out) :: af ! Final point
integer(i4b), intent(in) :: Nsteps, Nplt ! Number of intermediate points to be exported
interface
	subroutine SetNlin(a,N_a)
	use nrtype
	implicit none
	complex(dpc), dimension(:), intent(in) :: a
	complex(dpc), dimension(:), intent(out) :: N_a
	end subroutine
end interface
! Simple driver for etdrk4Diag. The functions f0-3,e,e2 are precomputed
! in calling routine.


integer(i4b) :: d, Npnt, plrt, i,j, dumm 
complex(dpc), dimension(size(ai)) :: a
real(dp) :: t,h

d=2*(assert_eq(size(f3), size(a),'etdrk4Diag-a')-1)

t=ti
a=ai

h=abs(tf-ti)/Nsteps

!print *,"t",Nsteps,Nsteps*h,tf
plrt=ceiling(Real(NSteps,dp)/Real(Nplt,dp)) ! Calculate after how many steps taken we should export values

if (allocated(tSt)) deallocate(tSt) !Clear out old stored variables if necessary.
if (allocated(aSt)) deallocate(aSt)
allocate(tSt(Nplt+1)) !Allocate storage for saved values.
allocate(aSt(Nplt+1,size(ai)))

j=1
tSt(j)=ti
aSt(1,:)=ai

do i=1,Nsteps
	call etdrk4Diag(a,h,a,f0,f1,f2,f3,e,e2,SetNlin)
	t=ti+i*h
	if (mod(i,plrt) == 0) then ! export some value
		j=j+1
		tSt(j)=t
		aSt(j,:)=a
!		print *,i,j,t,sum(abs(a))
	end if
end do

close(10)

af=a ! return final point if required

end subroutine
