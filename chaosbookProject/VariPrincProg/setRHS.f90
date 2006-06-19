SUBROUTINE setRHS(N,d,m1,m2,av,lambda,U,V,ch,b)

USE nrtype
USE interfcs

IMPLICIT NONE


INTEGER(I4B), INTENT(IN) :: N,d,m1,m2
REAL(DP), INTENT(IN)	:: av(:), lambda, U(:,:), V(:,:), ch
REAL(DP), INTENT(OUT)	::  b(:)

! Sets the rhs of the finite difference scheme. av contains the phase space
! coordinates in a column vector form (see main). 


REAL(DP) :: t, vel(d), q(size(b))
INTEGER(I4B) :: i

b=0
t=0 !!! if time is not explicit in derivs
vel=0


do i=1,N ! Form the right hand side vector, lambda*v part
	call derivs(t,av((i-1)*d+1:i*d),vel)
	b((i-1)*d+1:i*d)=lambda*vel
end do


! Form the right hand side vector, D.x part (without cyclic terms)
call banmul(FDM(N,d,m1,m2,ch),m1,m2,av,q) 




!Apply correction in q and form the RHS b.
b = ( b - ( q + MatMul(U,MatMul(Transpose(V),av)) ) )

END SUBROUTINE