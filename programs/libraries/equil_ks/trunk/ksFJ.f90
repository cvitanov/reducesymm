SUBROUTINE ksFJ(bc,fvec,fjac)
USE nrtype
use ifc_equil_ks
use nrutil, only:assert_eq

IMPLICIT NONE

include "fftw3.f"

real(dp), DIMENSION(:), INTENT(IN) :: bc
real(dp), DIMENSION(:), INTENT(OUT) :: fvec
real(dp), DIMENSION(:,:), INTENT(OUT) :: fjac
!!!!
complex(dpc), dimension(size(bc)/2) :: a
complex(dpc), dimension(size(bc)/2+1) :: adum 
complex(dpc), dimension(size(bc)/2+1) :: N_adum
complex(dpc), dimension(size(bc)/2) :: fvec_c !, fvecA
real(dpc), dimension(d/2,d/2):: jcc, jbb, jbc, jcb
integer(i4b):: ndum,k,j
real(dp), dimension(size(bc)) :: v 
integer(i8b) :: invplan, plan ! needed by fftw3
real(dp), dimension(size(bc)/2) :: q,lin
real(dp), dimension(size(fjac,1)):: wR,wI
complex(dpc), dimension(size(fjac,1)) :: W
complex(dpc), dimension(size(fjac,1),size(fjac,1)) :: VL, VR
integer(i4b) :: INFO

ndum=assert_eq(d,size(bc),size(fvec),'SetNlin1')
ndum=assert_eq(ndum,size(fjac,1),size(fjac,2),'SetNlin2')

! a does not include the a_0 coefficient
a=(0,0)
a=bc(1:d/2)+ ii*bc(d/2+1:d)
! 

adum=(0,0)
adum(2:size(adum))=a
call SetNlin_KS(adum,N_adum)

do k=1,d/2
	q(k)=k/L
	lin(k) = (1-(q(k))**2)*(q(k))**2
        fvec_c(k) = lin(k)*a(k) + N_adum(k+1)
end do


fvec(1:d/2)=real(fvec_c)
fvec(d/2+1:d)=aimag(fvec_c)

print *,"fvec", sum(abs(fvec))


jcc=0.0_dp
jbb=0.0_dp
jbc=0.0_dp
jcb=0.0_dp
fjac=0.0_dp

! Calculate Matrix of Variations(Jacobian)
!! calculate d\dot{c}/dc submatrix
do k=1,d/2
	jcc(k,k) = lin(k)
	do j=1,k-1
		jcc(k,j)=jcc(k,j)-2*q(k)*aimag(a(k-j))
	end do
	do j=k+1,d/2
		jcc(k,j)=jcc(k,j)+2*q(k)*aimag(a(j-k))
	end do
	do j=1,d/2-k
		jcc(k,j)=jcc(k,j)+2*q(k)*aimag(a(k+j))
	end do
end do
!! calculate d\dot{b}/db submatrix
do k=1,d/2
	jbb(k,k) = lin(k)
	do j=1,k-1
		jbb(k,j)=jbb(k,j)-2*q(k)*aimag(a(k-j))
	end do
	do j=k+1,d/2
		jbb(k,j)=jbb(k,j)+2*q(k)*aimag(a(j-k))
	end do
	do j=1,d/2-k
		jbb(k,j)=jbb(k,j)-2*q(k)*aimag(a(k+j))
	end do
end do
!! calculate d\dot{b}/dc submatrix
do k=1,d/2
	do j=1,k-1
		jbc(k,j)=-2*q(k)*real(a(k-j))
	end do
	do j=k+1,d/2
		jbc(k,j)=-2*q(k)*real(a(j-k))
	end do
	do j=1,d/2-k
		jbc(k,j)=jbc(k,j)+2*q(k)*real(a(k+j))
	end do
end do
!! calculate d\dot{c}/db submatrix
do k=1,d/2
	do j=1,k-1
		jcb(k,j)=2*q(k)*real(a(k-j))
	end do
	do j=k+1,d/2
		jcb(k,j)=jcb(k,j)+2*q(k)*real(a(j-k))
	end do
	do j=1,d/2-k
		jcb(k,j)=jcb(k,j)+2*q(k)*real(a(k+j))
	end do
end do

fjac(1:d/2,1:d/2)=jbb
fjac(1:d/2,d/2+1:d)=jbc
fjac(d/2+1:d,1:d/2)=jcb
fjac(d/2+1:d,d/2+1:d)=jcc

END SUBROUTINE