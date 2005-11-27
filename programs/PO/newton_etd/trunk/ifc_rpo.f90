module ifc_rpo

interface
	subroutine setLHM(f,J,v,q,R,diagk,LHM)
		use nrtype 
		implicit none	
		complex(dpc), dimension(:,:), intent(inout) :: J
		complex(dpc), dimension(:), intent(in) ::  f, v, q, R
		real(dp), dimension(:), intent(in) :: diagk
		complex(dpc), dimension(:,:), intent(out) :: LHM 
	end subroutine
end interface

interface
	subroutine setRHS(diff,RHS)
		use nrtype
		implicit none
		complex(dpc),  dimension(:), intent(in) :: diff
		complex(dpc),  dimension(:), intent(out) :: RHS
	end subroutine
end interface

interface
	subroutine newtonPOetdrk4(ai,T,kappa,q,tol,maxIter,h,diagk,f0,f1,f2,f3,e,e2,SetLin,SetNlin,SetANdiag,conv,J)
		use nrtype
		implicit none
		complex(dpc), intent(inout):: ai(:)
		real(dp), intent(inout) :: T		! Guess and calculated period
		real(dp), intent(inout) :: kappa	! Guess kappa in, computed kappa out		
		complex(dpc), intent(in) :: q(:) 	! Vector normal to Poincare hyperplane
		real(dp), intent(in) :: tol		
		integer(i4b), intent(in) :: maxIter	! Maximum number of iterations to be attempted
		real(dp), intent(in) :: h
		real(dp), dimension(:), intent(in) :: diagk	! Diag[ik/tilde{L}]	
		real(dp), dimension(:),intent(in) :: f0,f1,f2,f3,e,e2 ! precomputed functions of the linear operator
		integer(i4b), intent(out) :: conv	! If converged conv -> 1, if not conv -> 0
		complex(dpc), intent(out) :: J(:,:)
		interface
			subroutine SetLin(Lin)
				use nrtype
				implicit none
				real(dp), dimension(:), intent(out) :: Lin
			end subroutine
		end interface
		interface
			subroutine SetNlin(a,N_a)
			use nrtype
			implicit none
			complex(dpc), dimension(:), intent(in) :: a
			complex(dpc), dimension(:), intent(out) :: N_a
			end subroutine
		end interface
		interface
			subroutine SetAndiag(a,Andiag)
			use nrtype
			implicit none
			complex(dpc), dimension(:), intent(in) :: a
			complex(dpc), dimension(:,:), intent(out) :: Andiag
			end subroutine
		end interface
	end subroutine
end interface

end module