SUBROUTINE ksFJ(bc,fvec,fjac)
USE nrtype
use parameters
use nrutil, only:assert_eq

IMPLICIT NONE

include "fftw3.f"

real(dp), DIMENSION(:), INTENT(IN) :: bc
real(dp), DIMENSION(:), INTENT(OUT) :: fvec
real(dp), DIMENSION(:,:), INTENT(OUT) :: fjac
!!!!
complex(dpc), dimension(size(bc)/2) :: a
complex(dpc), dimension(size(bc)/2+1) :: adum 
complex(dpc), dimension(size(bc)/2) :: N_a
complex(dpc), dimension(size(bc)/2+1) :: N_adum
complex(dpc), dimension(size(bc)/2) :: fvec_c
real(dpc), dimension(d/2,d/2):: jcc, jbb, jbc, jcb
integer(i4b):: ndum,k ,j
real(dp), dimension(size(bc)) :: v 
integer(i8b) :: invplan, plan ! needed by fftw3
real(dp), dimension(size(bc)/2) :: q,lin

ndum=assert_eq(d,size(bc),size(fvec),'SetNlin1')
ndum=assert_eq(ndum,size(fjac,1),size(fjac,2),'SetNlin2')

! a does not include the a_0 coefficient
a=(0,0)
a=bc(1:size(bc)/2)+ ii*bc(size(bc)/2+1:size(bc))

adum=(0,0)
adum(2:size(a))=a
call dfftw_plan_dft_c2r_1d(invplan,d,adum,v,FFTW_ESTIMATE)
call dfftw_execute(invplan)
call dfftw_destroy_plan(invplan)
v=v**2
call dfftw_plan_dft_r2c_1d(plan,d,v,N_adum,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
N_a=N_adum(2:size(N_adum))/size(N_adum)

do k=1,d/2
	q(k)=k/L
	lin(k) = (1-(q(k))**2)*(q(k))**2
end do

do k=1,d/2 
	fvec_c(k) = lin(k)*a(k) + ii*q(k)*N_a(k) 
end do

fvec(1:d/2)=real(fvec_c)
fvec(d/2+1:d)=aimag(fvec_c)

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
		jcc(k,j)=jbb(k,j)+2*q(k)*aimag(a(k+j))
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
		jbc(k,j)=jbb(k,j)+2*q(k)*real(a(k+j))
	end do
end do
!! calculate d\dot{c}/db submatrix
do k=1,size(fvec)
	do j=1,k-1
		jcb(k,j)= 2*q(k)*real(a(k-j))
	end do
	do j=k+1,d/2
		jbc(k,j)= 0
	end do
end do

fjac(1:d/2,1:d/2)=jbb
fjac(1:d/2,d/2+1:d)=jbc
fjac(d/2+1:d,1:d/2)=jcb
fjac(d/2+1:d,d/2+1:d)=jcc



END SUBROUTINE