program kursiv

use nrtype
use parameters
use ifc_integr

implicit none

include "fftw3.f"


real(dp), dimension(d/2+1) :: Lin, f0,f1,f2,f3,e,e2
real(dp), dimension(d) :: u,v
integer(i4b) :: k,i
complex(dpc), dimension(d/2+1) :: ai,af
integer(i8b) :: plan,invplan
real(dp) :: rnd

interface
	subroutine SetLin_KS(Lin)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(out) :: Lin
	end subroutine
end interface
interface
	subroutine SetNlin_KS(a,N_a)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in) :: a
		complex(dpc), dimension(:), intent(out) :: N_a
	end subroutine
end interface


call SetLin_KS(Lin)

call etdrk4DiagPrefactors(Lin,h,R,M,f0,f1,f2,f3,e,e2)

u=0.0_dp
do i=1,d
	call Random_number(rnd)
	u(i)=(1+0.5*rnd)*twopi_d*real(i,dp)/real(d,dp)
end do

u=cos(u)*(1.0_dp+sin(u))

! open(10,file='inp.dat')
! read(10,*) u
! close(10)

call dfftw_plan_dft_r2c_1d(plan,d,u,ai,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)

!print *,ai(1),ai(2)

ai(1)=(0.0,0.0)

call etdrk4DiagDriverS(ti,ai,h,tf,af,f0,f1,f2,f3,e,e2,Nplt,SetNlin_KS)

!print *,size(aSt,1),size(aSt,2)

open(9,file='ksu.dat')
open(17,file='ksa.dat')
do i=1,size(aSt,1)
	write(17,*) real(aSt(i,2)),imag(aSt(i,2)),real(aSt(i,3))
!	print *,aSt(i,2)
        call dfftw_plan_dft_c2r_1d(invplan,d,aSt(i,:),v,FFTW_ESTIMATE)
	call dfftw_execute(invplan)
	call dfftw_destroy_plan(invplan)
!	print *,u
	v=v/size(v)
	write(9,frm_u) v
end do
close(9)
close(17)

!open(10,file='kst.dat')
!write(10,frm_t) tSt
!close(10)

open(11,file='inp.dat')
write(11,frm_u) v
close(11)


end program
