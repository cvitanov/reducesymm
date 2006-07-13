module ifc_po

interface
	subroutine setLHM(f,J,v,q,LHM)
		use nrtype 
		implicit none	
		complex(dpc), dimension(:,:), intent(in) :: J
		complex(dpc), dimension(:), intent(in) ::  f, v, q
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

interface newtonPOetdrk4
	subroutine newtonPOetdrk4_g(ai,T,q,tol,maxIter,h,f0,f1,f2,f3,e,e2,SetLin,SetNlin,SetANdiag,conv,J)
		use nrtype
		implicit none
		complex(dpc), intent(inout):: ai(:)
		real(dp), intent(inout) :: T		! Guess and calculated period
		complex(dpc), intent(in) :: q(:) 	! Vector normal to Poincare hyperplane
		real(dp), intent(in) :: tol		
		integer(i4b), intent(in) :: maxIter	! Maximum number of iterations to be attempted
		real(dp), intent(in) :: h
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
	end subroutine newtonPOetdrk4_g
	subroutine newtonPOetdrk4_c0(ai,T,q,tol,maxIter,h,f0,f1,f2,f3,e,e2,SetLin,SetNlin,SetANdiag,conv,J,c0)
		use nrtype
		implicit none
		complex(dpc), intent(inout):: ai(:)
		real(dp), intent(inout) :: T		! Guess and calculated period
		complex(dpc), intent(in) :: q(:) 	! Vector normal to Poincare hyperplane
		real(dp), intent(in) :: tol		
		integer(i4b), intent(in) :: maxIter	! Maximum number of iterations to be attempted
		real(dp), intent(in) :: h
		real(dp), dimension(:),intent(in) :: f0,f1,f2,f3,e,e2 ! precomputed functions of the linear operator
		integer(i4b), intent(out) :: conv	! If converged conv -> 1, if not conv -> 0
		complex(dpc), intent(out) :: J(:,:)
		integer(i4b), intent(in), optional:: c0 !
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
	end subroutine newtonPOetdrk4_c0
end interface

end module