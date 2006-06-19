FUNCTION FDM(N,d,m1,m2,ch)

USE nrtype

IMPLICIT NONE

INTEGER(I4B), INTENT(IN):: N,m1,m2,d
REAL(DP), INTENT(IN) :: ch
REAL(DP) :: FDM(N*d+1,m1+1+m2)

! Computes the matrix D, the discrete differentiation matrix.



INTEGER(I4B) :: i

FDM=0

! Constract the matrice of finite differences
do i=1,(N-1)*d !Write the 8I Matrix in the superdiagonal
	FDM(i,(m1+1)+d)=8
end do

do i=1,(N-2)*d !Write the -I Matrix in the superdiagonal
	FDM(i,(m1+1)+2*d)=-1
end do

do i=d+1,N*d !Write the -8I Matrix in the subdiagonal
	FDM(i,(m1+1)-d)=-8
end do

do i=2*d+1,N*d !Write the I Matrix in the subdiagonal
	FDM(i,1)=1
end do

FDM=ch*FDM ! Include the prefactor


END FUNCTION