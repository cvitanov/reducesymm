function DiagMul_rc(v,M,N)

use nrtype

implicit none

real(dp), dimension(N), intent(in) :: v
complex(dpc), dimension(N,N), intent(in) :: M
integer(i4b), intent(in) :: N
complex(dpc), dimension(N,N) :: DiagMul_rc
! Performs matrix multiplication of a diagonal NxN real matrix whose
! diagonal is stored in vector v with matrix M. Thus the result
! is diag(v).M .
! v: real
! M: complex

integer(i4b) :: i

do i=1,size(v)
	DiagMul_rc(i,:) = v(i)*M(i,:)
end do

end function DiagMul_rc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function DiagMul_cc(v,M,N)

use nrtype

implicit none

complex(dpc), dimension(N), intent(in) :: v
complex(dpc), dimension(N,N), intent(in) :: M
integer(i4b), intent(in) :: N
complex(dpc), dimension(N,N) :: DiagMul_cc
! Performs matrix multiplication of a diagonal NxN real matrix whose
! diagonal is stored in vector v with matrix M. Thus the result
! is diag(v).M .
! v: complex
! M: complex

integer(i4b) :: i



do i=1,size(v)
	DiagMul_cc(i,:) = v(i)*M(i,:)
end do

end function DiagMul_cc