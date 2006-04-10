SUBROUTINE ksFJ(a,fvec,fjac)
USE nrtype
use parameters
use nrutil, only:assert_eq

IMPLICIT NONE

include "fftw3.f"

complex(dpc), DIMENSION(:), INTENT(IN) :: a
complex(dpc), DIMENSION(:), INTENT(OUT) :: fvec
complex(dpc), DIMENSION(:,:), INTENT(OUT) :: fjac
!!!!
complex(dpc), dimension(size(a)) :: N_a
complex(dpc), dimension(size(a)+1) :: N_adum
integer(i4b):: ndum,k ,j
real(dp), dimension(2*(size(a))) :: v 
complex(dpc), dimension(size(a)+1) :: adum 
integer(i8b) :: invplan, plan ! needed by fftw3
real(dp), dimension(size(fvec)) :: q,lin

ndum=assert_eq(d,2*size(a),'SetNlin')
ndum=assert_eq(size(N_a),size(fvec),'SetNlin')

!print *,"v",size(v),"a",size(a),"N_adum", size(N_adum), "a_dum",size(adum),"fvec", size(fvec), "fjac", size(fjac), "N_a", size(N_a)

adum(2:size(adum))=a
call dfftw_plan_dft_c2r_1d(invplan,d,adum,v,FFTW_ESTIMATE)
call dfftw_execute(invplan)
call dfftw_destroy_plan(invplan)
v=v**2
call dfftw_plan_dft_r2c_1d(plan,d,v,N_adum,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
N_a=N_adum(2:size(N_adum))/size(N_adum)


do k=1,size(fvec)
	q(k)=k/L
	lin(k) = (1-(q(k))**2)*(q(k))**2
	print *,k,lin(k) 
end do

do k=1,size(fvec) 
	fvec(k) = lin(k)*a(k) + ii*q(k)*N_a(k) 
end do

! Calculate Matrix of Variations(Jacobian)
do k=1,size(fvec)
	fjac(k,k) = lin(k)
	do j=1,k-1
		fjac(k,j)=fjac(k,j)+2*ii*q(k)*a(k-j)
		if ( k == j ) print *,k,2*ii*q(k)*a(k-j)
	end do
	do j=k+1,size(fvec)
		fjac(k,j)=fjac(k,j)-2*ii*q(k)*a(j-k)
		if ( k == j ) print *,k,-2*ii*q(k)*a(j-k)
	end do
	do j=1,size(fvec)-k
		fjac(k,j)=fjac(k,j)-2*ii*q(k)*a(k+j)
		if ( k == j ) print *,k,-2*ii*q(k)*a(k+j)
	end do
	print *,"diag",k,fjac(k,k)
end do

open(20,file="rA.dat")
open(21,file="iA.dat")
	do k=1,size(fvec)
		do j=1,size(fvec)
			write(20,*) real(fjac(k,j))
			write(21,*) aimag(fjac(k,j))
		end do 
	end do
close(20)
close(21)

open(9,file='a0r.dat')
open(10,file='a0i.dat')
	write (9,*) real(fvec)
	write (10,*) aimag(fvec)
close(9)
close(10)


END SUBROUTINE