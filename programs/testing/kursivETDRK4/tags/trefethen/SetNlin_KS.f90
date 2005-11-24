subroutine SetNlin_KS(a,N_a)

use nrtype
use parameters
use nrutil, only:assert_eq


implicit none

include "fftw3.f"

complex(dpc), dimension(:), intent(in) :: a
complex(dpc), dimension(:), intent(out) :: N_a
! Returns the nonlinear part N_a in the KSe. Notice it is
! equal to the nonlinear operator acting on a. 
integer(i4b):: ndum,k 
real(dp), dimension(2*(size(a)-1)) :: v 
complex(dpc), dimension(size(a)) :: adum 
integer(i8b) :: invplan, plan ! needed by fftw3


ndum=assert_eq(d,2*(size(a)-1),2*(size(N_a)-1),'SetNlin')

adum=a
call dfftw_plan_dft_c2r_1d(invplan,d,adum,v,FFTW_ESTIMATE)
call dfftw_execute(invplan)
call dfftw_destroy_plan(invplan)
v=v**2
call dfftw_plan_dft_r2c_1d(plan,d,v,N_a,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
N_a=N_a/size(v)
do k=1,size(N_a) 
	N_a(k)=-0.5_dp*ii*(k-1)*N_a(k)/L
end do

end subroutine

