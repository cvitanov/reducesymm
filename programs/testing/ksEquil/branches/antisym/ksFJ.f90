SUBROUTINE ksFJ(aim,fvec,fjac)
USE nrtype
use parameters
use nrutil, only:assert_eq

IMPLICIT NONE

include "fftw3.f"

real(dp), DIMENSION(:), INTENT(IN) :: aim
real(dp), DIMENSION(:), INTENT(OUT) :: fvec
real(dp), DIMENSION(:,:), INTENT(OUT) :: fjac
interface
	subroutine SetNlin_KS(a,N_a)
		use nrtype
		implicit none
		complex(dpc), dimension(:), intent(in) :: a
		complex(dpc), dimension(:), intent(out) :: N_a
	end subroutine
end interface
!!!!
complex(dpc), dimension(size(aim)) :: a
complex(dpc), dimension(size(aim)+1) :: adum 
complex(dpc), dimension(size(aim)) :: N_a
complex(dpc), dimension(size(aim)+1) :: N_adum,N_af
complex(dpc), dimension(size(aim)) :: fvec_c
integer(i4b):: ndum,k ,j,m
real(dp), dimension(d) :: v 
integer(i8b) :: invplan, plan ! needed by fftw3
real(dp), dimension(d/2) :: q,lin

ndum=assert_eq(d,size(aim)*2,size(fvec)*2,'SetNlin1')
ndum=assert_eq(ndum,size(fjac,1)*2,size(fjac,2)*2,'SetNlin2')

! a does not include the a_0 coefficient
a=(0,0)
a=ii*aim

adum=(0,0)
adum(2:size(adum))=ii*aim

call SetNlin_KS(adum,N_adum)

do k=1,d/2
	q(k)=k/L
	lin(k) = (1-(q(k))**2)*(q(k))**2
end do

do k=1,d/2 
	fvec_c(k) =  lin(k)*ii*aim(k) + N_adum(k+1) 
end do

fvec=aimag(fvec_c)

print *,"fvec",sum(abs(fvec))

fjac=0.0_dp

! Calculate Matrix of Variations(Jacobian)
!! calculate d\dot{c}/dc submatrix
do k=1,d/2
	fjac(k,k) = lin(k)
	do j=1,k-1
		fjac(k,j)=fjac(k,j)-2*q(k)*aim(k-j)
	end do
	do j=k+1,d/2
		fjac(k,j)=fjac(k,j)+2*q(k)*aim(j-k)
	end do
	do j=1,d/2-k
		fjac(k,j)=fjac(k,j)+2*q(k)*aim(k+j)
	end do
end do



END SUBROUTINE