PROGRAM integr

USE nrtype
USE ifc_integr
USE ifc_util, ONLY: UnitMatrix
!USE F95_LAPACK, ONLY: LA_GEEV
USE parameters

IMPLICIT NONE

complex(dpc), DIMENSION(d,d) ::  J
REAL(DP) :: ti=0.0_dp, tf, h
INTEGER(I4B) :: nsteps=10000, Nplt=1000, i, M=32
complex(dpc), DIMENSION(d) :: yi
real(dp), dimension(d) :: yRi
!real(dp), allocatable :: WI(:), WR(:)
real(dp), dimension(d) :: Lin, f0,f1,f2,f3,e,e2
REAL(dp) :: R=1.0_dp


interface
	subroutine SetLin_roessler(Lin)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(out) :: Lin
	end subroutine
end interface
interface
	subroutine SetNlin_roessler(a,N_a)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in) :: a
		complex(dpc), dimension(:), intent(out) :: N_a
	end subroutine
end interface
interface
	subroutine SetANdiag_roessler(a,ANdiag)
	use nrtype
	implicit none
	complex(dpc), dimension(:), intent(in) :: a
	complex(dpc), dimension(:,:), intent(out) :: ANdiag
	end subroutine
end interface


INTEGER(I4B) :: idum, p=0

open(9,file='initial.dat')

read(9,*)  yi(1:d)

close(9)

open(17,file='period.dat')

read(17,*) tf

close(17)

h=(tf-ti)/nsteps


yi(d+1)=ti

p=0

print *,ti,tf
print *,yi

J=UnitMatrix(d)

call SetLin_Roessler(Lin)
call etdrk4DiagPrefactors(Lin,h,R,M,f0,f1,f2,f3,e,e2)

call etdrk4DiagJDriverS(ti,yi,J,h,tf,yi,J,f0,f1,f2,f3,e,e2,Nplt,SetNlin_roessler,SetANdiag_roessler)

!allocate(WR(size(J,1)),WI(size(J,1)))

!call LA_GEEV( J, WR, WI)

!print *, WR
!print *, WI



open(7,file="cycle.dat")
do i=1,size(aSt,1)
	write(7,format_label) Real(aSt(i,1:d))
end do
close(7)

open(9,file='J.dat')

do i=1,d
	write(9,format_label) Real(J(i,:))
end do

close(9)

!CLOSE(11)

!CLOSE(10) 


END PROGRAM