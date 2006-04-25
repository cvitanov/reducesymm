program POshootingTest

use nrtype
use parameters
use ifc_po
use ifc_integr
!use la_precision, only: wp => dp
!use f95_lapack
use ifc_util

implicit none

include "fftw3.f"


real(dp), dimension(d/2+1) :: Lin, f0,f1,f2,f3,e,e2
real(dp), dimension(d) :: u,v
integer(i4b) :: k,i, conv=0, c0
complex(dpc), dimension(d/2+1) :: ai,af,adum
complex(dpc), dimension(d+2) :: q
integer(i8b) :: plan,invplan
real(dp) ::  T
complex(dpc), dimension(d/2,d/2) ::  J

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
interface
	subroutine SetAndiag_KS(a,Andiag)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in) :: a
		complex(dpc), dimension(:,:), intent(out) :: Andiag
	end subroutine
end interface

call SetLin_KS(Lin)

call etdrk4DiagPrefactors(Lin,h,R,M,f0,f1,f2,f3,e,e2)

u=0.0_dp

open(10,file='guess.dat')
read(10,*) u
close(10)

open(11,file='Tguess.dat')
read(11,*) T
close(11)

q=(0.0,0.0)

q(2)=(1.0,0.0)
q(d/2+1+2+2)=-(1.0,0.0)

J=UnitMatrix(d/2+1)

call dfftw_plan_dft_r2c_1d(plan,d,u,ai,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)

ai=ai/size(u)

call newtonPOetdrk4(ai,T,q,tol,maxIter,h,f0,f1,f2,f3,e,e2,SetLin_KS,SetNlin_KS,SetANdiag_KS,conv,J,c0=1)

open(9,file='ksu.dat')
open(12,file='ksa.dat')
do i=1,size(aSt,1)
!	print *,aSt(i,2)
	write(12,'(4F15.10)') real(aSt(i,3)),imag(aSt(i,3)),real(aSt(i,5)),imag(aSt(i,5))
	adum=aSt(i,:)
        call dfftw_plan_dft_c2r_1d(invplan,d,adum,v,FFTW_ESTIMATE)
	call dfftw_execute(invplan)
	call dfftw_destroy_plan(invplan)
!	print *,u
!	v=v/size(v)
	write(9,frm_u) v
end do
close(9)

!open(10,file='kst.dat')
!write(10,frm_t) tSt
!close(10)

open(11,file='inp.dat')
write(11,frm_u) v
close(11)


end program
