program kursiv

use nrtype
use parameters
use ifc_integr

implicit none

include "fftw3.f"


real(dp), dimension(d/2+1) :: Lin, f0,f1,f2,f3,e,e2
real(dp), dimension(d) :: u,v
integer(i4b) :: k,i
complex(dpc), dimension(d/2+1) :: ai,af,adum
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

open(10,file='guess.dat')
read(10,*) u
close(10)

call dfftw_plan_dft_r2c_1d(plan,d,u,ai,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)

!print *,ai(1),ai(2)

ai(1)=(0.0,0.0)

do i=int(size(ai)/2),size(ai)
   ai(i)=(0.0,0.0)
end do

ai=ai/size(u)

call etdrk4DiagDriverS(ti,ai,h,tf,af,f0,f1,f2,f3,e,e2,Nplt,SetNlin_KS)

!print *,size(aSt,1),size(aSt,2)

open(9,file='ksu.dat')
open(17,file='ksa.dat')
open(18,file='ksaAll.dat')
do i=1,size(aSt,1)
	write(17,'(4F15.10)') real(aSt(i,3)),imag(aSt(i,3)),real(aSt(i,5)),imag(aSt(i,5))
!	print *,aSt(i,2)
	adum=aSt(i,:)
        call dfftw_plan_dft_c2r_1d(invplan,d,adum,v,FFTW_ESTIMATE)
	call dfftw_execute(invplan)
	call dfftw_destroy_plan(invplan)
!	print *,v
!	v=v/size(v)
	write(9,frm_u) v
end do
close(9)
close(17)
close(18)

!open(10,file='kst.dat')
!write(10,frm_t) tSt
!close(10)

open(11,file='inp.dat')
write(11,frm_u) v
close(11)


end program
