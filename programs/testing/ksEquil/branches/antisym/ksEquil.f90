Program ksEquil

use nrtype
use parameters
use ifc_newt

implicit none

include "fftw3.f"

real(dp), dimension(d) :: v
complex(dpc), dimension(d/2+1) :: a
real(dp), dimension(d/2) :: aim
integer(i8b) :: invplan, plan ! needed by fftw3
integer(i4b) :: k,i
!!!!
interface
	subroutine ksFJ(bc,fvec,fjac)
		USE nrtype
		implicit none
		real(dp), DIMENSION(:), INTENT(IN) :: bc
		real(dp), DIMENSION(:), INTENT(OUT) :: fvec
		real(dp), DIMENSION(:,:), INTENT(OUT) :: fjac
	end subroutine
end interface

!!!!


!open(19,file='guess.dat')

!read(19,*) v

!close(19)

do i=1,d
	!print *,2*PI_D*Real(i)/Real(d),0.1*sin(2*PI_D*Real(i)/Real(d))
	v(i)=0.1*sin(2*PI*Real(i-1)/Real(d)) !+ 0.2*cos(2*PI*i/d) 
end do

!print *,v

call dfftw_plan_dft_r2c_1d(plan,d,v,a,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
a=a/size(v)

open(23,file='a0.dat')

do k=1,size(a)
	write(23,*) a(k)
end do

aim=aimag(a(2:size(a)))

print *,aim(1)

call mnewt(ntrial,aim,tola,tolf,ksFJ)

print *, "converged?"

open(23,file='equil.dat')

do k=1,size(aim)
	write(23,*) aim(k)
end do

close(23)


end program