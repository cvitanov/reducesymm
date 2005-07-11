subroutine etdrk4DiagPrefactors(Lin,h,R,M,f0,f1,f2,f3,e,e2)

use nrtype
use nrutil, only: assert_eq

implicit none

real(dp), dimension(:), intent(in) :: Lin ! The linear operator
real(dp), intent(in) :: h, R ! stepsize, radious of contour around each eigenvalue
integer(i4b),intent(in) :: M ! number of points to use in integration
real(dp), dimension(:),intent(out) :: f0,f1,f2,f3,e,e2 ! Various functions of Lin
! Precompute the functions of linear operator Lin required by etdrk4Diag.
! For the diagonal case the integration reduces to calculating the mean value of the
! function on the points of the contour for each eigenvalue of the
! linear operator. Lin is thus represented by a vector containing each
! eigenvalues. 
integer(i4b) :: ndum, i, k
complex(dpc), dimension(M) :: rootsUnity, fdum


ndum=assert_eq(size(Lin),size(e),size(e2),'etdrk4PrefactorsDiag-e')
ndum=assert_eq(ndum,size(f0),size(f1),'etdrk4PrefactorsDiag-0-1')
ndum=assert_eq(ndum,size(f2),size(f3),'etdrk4PrefactorsDiag-2-3')

do i=1,M
	rootsUnity(i)=exp(ii*PI_D*(i-0.5_dp)/M)
end do

e=exp(h*Lin) ! No problem with these two guys, calculate directly.
e2=exp(h*Lin/2)

do k=1,size(Lin) ! For each eigenvalue
	do i=1,M ! add contribution of each point on the contour
		fdum(i) = h*Lin(k)+R*rootsUnity(i)
	end do
	!and finally take mean values for each scalar function fi, to find fi(k) for eigenvalue Lin(k)
	f0(k) = h*Real(sum( (exp(fdum/2.0_dp)-1)/fdum )/size(fdum) )
	f1(k) = h*Real(sum( (-4-fdum+exp(fdum)*(4-3*fdum+fdum**2) )/fdum**3 )/size(fdum) )
	f2(k) = h*Real(sum( 2*(2+fdum+exp(fdum)*(-2+fdum) )/fdum**3 )/size(fdum) )
	f3(k) = h*Real(sum( (-4-3*fdum-fdum**2+exp(fdum)*(4-fdum) )/fdum**3 )/size(fdum) )
end do

end subroutine