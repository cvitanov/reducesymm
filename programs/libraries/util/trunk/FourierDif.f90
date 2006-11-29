subroutine FourierDif(v,vx,L,n)

use nrtype
use nrutil

implicit none

include "fftw3.f"

real(dp), dimension(:), intent(in) :: v
real(dp), dimension(:), intent(out) :: vx
real(dp), intent(in) :: L
integer(i4b), intent(in) :: n
! Given a tabulated set of values v of a function v(x) periodic in [0,2*pi*L],
! The subroutine calculates the numerical derivative of order n
! vx = d^n v/dx^n by transforming to Fourier space.
integer(i8b) :: invplan, plan ! needed by fftw3
complex(dp), dimension(:), allocatable :: a
integer(i4b) :: d,k

d=assert_eq(size(v),size(vx),'FourierDif')

allocate(a(d/2+1))

call dfftw_plan_dft_r2c_1d(plan,d,v,a,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
a=a/size(v)

do k=1,d/2+1
	a(k)=((real(k-1,dp)/L)**n)*(ii**n)*a(k)
end do

call dfftw_plan_dft_c2r_1d(invplan,d,a,vx,FFTW_ESTIMATE)
call dfftw_execute(invplan)
call dfftw_destroy_plan(invplan)



end subroutine