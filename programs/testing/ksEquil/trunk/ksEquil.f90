Program ksEquil

use nrtype
use parameters
use ifc_newt

implicit none

include "fftw3.f"

real(dp), dimension(d) :: v
complex(dpc), dimension(d/2+1) :: a
real(dp), dimension(d) :: bc
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


open(19,file='guess2.dat')
 
	read(19,*) v
 
close(19)

! do i=1,d
! !	print *,2*PI_D*Real(i)/Real(d),0.1*sin(2*PI_D*Real(i)/Real(d))
! 	v(i)=0.2*sin(2*PI*Real(i-1)/Real(d)) !+ 0.2*cos(2*PI*(i-1)/d) 
! end do

!print *,v

call dfftw_plan_dft_r2c_1d(plan,d,v,a,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
a=a/size(v)

open(23,file='a0.dat')

do k=1,size(a)
	write(23,*) a(k)
end do

print *, "a",size(a), "bc", size(bc)

bc(1:d/2)=real(a(2:size(a)))
bc(d/2+1:d)= aimag(a(2:size(a)))

call mnewt(ntrial,bc,tolbc,tolf,ksFJ)

print *, "converged?"

open(23,file='equil.dat')

a=(0,0)

a(2:size(a))=bc(1:d/2)+ii*bc(d/2+1:d)

do k=1,size(a)
	write(23,*) a(k)
end do

close(23)

call dfftw_plan_dft_c2r_1d(invplan,d,a,v,FFTW_ESTIMATE)
call dfftw_execute(invplan)
call dfftw_destroy_plan(invplan)

open(24,file='Uequil.dat')

do k=1,size(v)
	write(24,*) v(k)
end do

close(24)

end program