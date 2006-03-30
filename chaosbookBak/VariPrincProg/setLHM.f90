	SUBROUTINE setLHM(N,d,m1,m2,av,lambda,ch,LHM)
	
	USE interfcs
	USE nrtype

	IMPLICIT NONE

	INTEGER(I4B), INTENT(IN) :: N,d,m1,m2
	REAL(DP), INTENT(IN) :: av(:), lambda, ch 
	REAL(DP), INTENT(OUT) :: LHM(:,:)
	
	! Creates the LHM in the finite difference scheme. av contains the phase space
	! coordinates in a column vector form (see main).

	INTEGER(I4B) :: i,m,j
	REAL(DP) :: MVdum(d,d),t,vel(d)

	t=0 !!! if derivs does not contain time explicity
	LHM=0

	! Assign values to the left hand side matrix.
	LHM=FDM(N,d,m1,m2,ch)	! Copy the values from the finite difference matrix

	LHM(N*d+1,(m1+1)-d)=1 ! Fix the N'th point of the Loop to be on the Poincare section a1=0. 

	do i=1,d !Write the velocity elements in the superdiagonal
		call derivs(t,av((N-1)*d+1:N*d),vel)
		LHM( (N*d+1)-i,(m1+1)+i)=-vel(d+1-i)
	end do
	do i=1,d
		call derivs(t,av((N-2)*d+1:(N-1)*d),vel)
		LHM( ((N-1)*d+1)-i,(m1+1)+d+i)=-vel(d+1-i)
	end do

	do m=1,d !Write the diagonal and superdiagonal elements involving the matrix of variations
		do i=1,N
			MVdum=MatVar(av((i-1)*d+1:i*d),d)
			do j=1,d-(m-1)
				LHM((i-1)*d+j,(m1+1)+(m-1)) = -lambda*MVdum(j,j+(m-1))
			end do
		end do
	end do


	do m=2,d !Write the subdiagonal elements involving the matrix of variations
		do i=1,N
			MVdum=MatVar(av((i-1)*d+1:i*d),d)
			do j=1,d-(m-1)
				LHM((i-1)*d+j+(m-1),(m1+1)-(m-1)) = -lambda*MVdum(j+(m-1),j)
			end do
		end do
	end do

	END SUBROUTINE



