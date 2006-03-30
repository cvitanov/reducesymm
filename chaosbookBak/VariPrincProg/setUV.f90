SUBROUTINE setUV(N,d,ch,U,V)

USE nrtype
USE interfcs
IMPLICIT NONE

INTEGER(I4B), INTENT(IN) :: N,d
REAL(DP), INTENT(IN) :: ch
REAL(DP), INTENT(OUT) :: U(:,:),V(:,:)


! Form the correction matrices to include cyclic terms, but not the velocity
! boundary terms. 

INTEGER(I4B) i

U=0
V=0

U((N-2)*d+1:(N-1)*d,1:d)=-UnitMatrix(d)
U((N-1)*d+1:N*d,1:d)=8*UnitMatrix(d)
U((N-1)*d+1:N*d,d+1:2*d)=-UnitMatrix(d)
U(1:d,2*d+1:3*d)=UnitMatrix(d)
U(1:d,3*d+1:4*d)=-8*UnitMatrix(d)
U(d+1:2*d,3*d+1:4*d)=UnitMatrix(d)

U=ch*U

do i=1,2*d
	V(i,i)= 1 
	V((N-2)*d+i,2*d+i)=1
end do

END subroutine