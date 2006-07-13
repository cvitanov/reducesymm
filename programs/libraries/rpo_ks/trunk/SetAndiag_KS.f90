subroutine SetANdiag_KS(a,Andiag)

use nrtype
use ifc_rpo_ks
use nrutil, only:assert_eq

implicit none

complex(dpc), dimension(:), intent(in) :: a
real(dp), dimension(:,:), intent(out) :: Andiag
!
!
real(dp), dimension(size(a),size(a)):: Acc, Abb, Abc, Acb
integer(i4b):: d,k,j
real(dp), dimension(size(a)) :: q

d=assert_eq(2*size(a),size(Andiag,1),size(Andiag,2),'SetANdiag')

Acc=0.0_dp
Abb=0.0_dp
Abc=0.0_dp
Acb=0.0_dp

do k=1,d/2
	q(k)=k/L
end do

! Calculate Matrix of Variations(Jacobian) - Non diagonal Part ONLY!
!! calculate d\dot{c}/dc submatrix
do k=1,d/2
	do j=1,k-1
		Acc(k,j)=Acc(k,j)-2.0_dp*q(k)*aimag(a(k-j))
	end do
	do j=k+1,d/2
		Acc(k,j)=Acc(k,j)+2.0_dp*q(k)*aimag(a(j-k))
	end do
	do j=1,d/2-k
		Acc(k,j)=Acc(k,j)+2.0_dp*q(k)*aimag(a(k+j))
	end do
end do
!! calculate d\dot{b}/db submatrix
do k=1,d/2
	do j=1,k-1
		Abb(k,j)=Abb(k,j)-2.0_dp*q(k)*aimag(a(k-j))
	end do
	do j=k+1,d/2
		Abb(k,j)=Abb(k,j)+2.0_dp*q(k)*aimag(a(j-k))
	end do
	do j=1,d/2-k
		Abb(k,j)=Abb(k,j)-2.0_dp*q(k)*aimag(a(k+j))
	end do
end do
!! calculate d\dot{b}/dc submatrix
do k=1,d/2
	do j=1,k-1
		Abc(k,j)=Abc(k,j)-2.0_dp*q(k)*real(a(k-j))
	end do
	do j=k+1,d/2
		Abc(k,j)=Abc(k,j)-2.0_dp*q(k)*real(a(j-k))
	end do
	do j=1,d/2-k
		Abc(k,j)=Abc(k,j)+2.0_dp*q(k)*real(a(k+j))
	end do
end do
!! calculate d\dot{c}/db submatrix
do k=1,d/2
	do j=1,k-1
		Acb(k,j)=Acb(k,j)+2.0_dp*q(k)*real(a(k-j))
	end do
	do j=k+1,d/2
		Acb(k,j)=Acb(k,j)+2.0_dp*q(k)*real(a(j-k))
	end do
	do j=1,d/2-k
		Acb(k,j)=Acb(k,j)+2.0_dp*q(k)*real(a(k+j))
	end do
end do

Andiag(1:d/2,1:d/2)=Abb
Andiag(1:d/2,d/2+1:d)=Abc
Andiag(d/2+1:d,1:d/2)=Acb
Andiag(d/2+1:d,d/2+1:d)=Acc

end subroutine
