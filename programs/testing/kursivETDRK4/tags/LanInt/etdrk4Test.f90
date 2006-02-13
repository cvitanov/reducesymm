program kursiv

use nrtype
use parameters
use ifc_integr

implicit none

include "fftw3.f"


real(dp), dimension(d/2+1) :: Lin, f0,f1,f2,f3,e,e2
real(dp), dimension(d) :: u,v
integer(i4b) :: k,i, skip=10, j
complex(dpc), dimension(d/2+1) :: ai,af
integer(i8b) :: plan,invplan
!real(dp) :: sgn

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
! do i=1,d
! 	u(i)=twopi_d*real(i,dp)/real(d,dp)
! end do
! 
! u=cos(u)*(1.0_dp+sin(u))

!open(10,file='inp.dat')
!read(10,*) u
!close(10)

!call dfftw_plan_dft_r2c_1d(plan,d,u,ai,FFTW_ESTIMATE)
!call dfftw_execute(plan)
!call dfftw_destroy_plan(plan)

!print *,ai(1),ai(2)

! ai=(0.0_dp,0.0_dp)
! ai(2)=(0.0_dp,0.06)
! ai(4)=(0.0_dp,0.075_dp) 
! ai(5)=(0.0_dp,-0.11_dp)
! ai(7)=(0.0_dp,-0.05_dp)
! ai(8)=(0.0_dp,-0.029_dp)

!do j=1,Niter

open(9,file='a0.dat')
read(9,*) ai
close(9)


open(8,file='ksF.dat')

do j=1,Niter
	print *,j
	call etdrk4DiagDriverS_a(ti,ai,h,tf,af,f0,f1,f2,f3,e,e2,Nplt,SetNlin_KS)
	do i=1,size(aSt,1)
		if ((aimag(aSt(i,2)) < 0.06_dp) .and. (aimag(aSt(i+1,2)) >=  0.06_dp) ) then
			write(8,'(3F15.5)') aimag(aSt(i,2)), aimag(aSt(i,4)),aimag(aSt(i,7))  
		!	print *, aSt(i,1)
		end if
	end do
	ai=af
end do

close(8)

!open(14,file='a0.dat')
!write(14,frm_a) aimag(aSt(size(aSt,1),:))
!close(14)

! open(9,file='ksu.dat')
! do i=1,size(aSt,1),skip
! !	print *,aSt(i,2)
!         call dfftw_plan_dft_c2r_1d(invplan,d,aSt(i,:),v,FFTW_ESTIMATE)
! 	call dfftw_execute(invplan)
! 	call dfftw_destroy_plan(invplan)
! !	print *,u
! 	v=v/size(v)
! 	write(9,frm_u) v
! end do
! close(9)

!open(10,file='kst.dat')
!write(10,frm_t) tSt
!close(10)

!open(11,file='inp.dat')
!write(11,frm_u) v(size(v))
!close(11)


end program
