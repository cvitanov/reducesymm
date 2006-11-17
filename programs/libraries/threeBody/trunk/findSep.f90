subroutine findSep_TBL(E,Z,xo,xmax,pUp,pDown,tol,Npoints,MaxAttempts,eps,h1,hmin,ti,taui,tauf,EOM_TBL)

use nrtype
use ifc_threeBody

implicit none

real(dp), intent(in) :: E, Z, xo,xmax
real(dp), dimension(:), intent(inout) ::  pUp, pDown
real(dp), intent(in) :: tol 
integer(i4b), intent(in) :: Npoints, MaxAttempts
real(dp), intent(in) :: eps, h1, hmin
real(dp), intent(in) :: ti, taui, tauf
interface 
	subroutine EOM_TBL(tau,f,v)
		use nrtype
		implicit none
		real(dp), dimension(:), intent(in) :: f
		real(dp), intent(in) :: tau
		real(dp), dimension(:), intent(out) :: v
	end subroutine
end interface

!!
real(dp), dimension(2) :: xi
real(dp), dimension(7) :: fi
real(dp) :: diff
real(dp), dimension(2) :: Q,P,r
real(dp) :: Et,hyperR
logical(lgt) :: pChangeUp = .false., pChangeDown = .false.
integer(i4b) :: i,j=0
real(dp) :: dx, incr


xi(1)=xo

dx= (xmax-xo)/Npoints

open(19,file='manifold.dat')
open(20,file='check.dat')

do i=1, Npoints
	print *,"Point #",i,xi(1)
	diff=pUp(1)-pDown(1)
	do while( diff > tol)
		!!! First integrate large momentum ic.
		fi=0.0_dp
		call findIC(xi,pUp,E,Z)
		call CartToKSMcG_TBL(xi,pUp,E,Q,P,r,Et,hyperR)
		call vectorize(Q,P,Et,hyperR,ti,fi) 
		call TBLint_sep(fi,taui,tauf,eps,h1,hmin,EOM_TBL,TBrkqs,pChangeUp) 
		!!! Then integrate small momentum ic.
		fi=0.0_dp
		call findIC(xi,pDown,E,Z)
		call CartToKSMcG_TBL(xi,pDown,E,Q,P,r,Et,hyperR)
		call vectorize(Q,P,Et,hyperR,ti,fi)
		call TBLint_sep(fi,taui,tauf,eps,h1,hmin,EOM_TBL,TBrkqs,pChangeDown)
		print *,"diff",diff,"pUp",pUp,"pDown",pDown
		if ( pChangeUp == .true. ) then !! Bad guess of pUp(1)
			incr =  2.0_dp*(pUp(1)-pDown(1))
			pDown(1) = pUp(1) 
			pUp(1)=pUp(1)+incr
		else 	
			if (pChangeDown == .false.) then !! Bad guess of pDown(1)
				pUp(1)=pDown(1)
				pDown(1) = 0.5_dp*pDown(1)
			else	!  both correct
				pUp(1) = 0.5_dp*(pUp(1)+pDown(1))	
			end if		
		end if
		diff=pUp(1)-pDown(1)
		j=j+1
		if (j>maxAttempts) then
			print *, "Convergence to",diff," after", MaxAttempts,"iterations."
			stop
		end if
	end do
	write(19,'(2F15.10)') xi(1),0.5_dp*(pUp(1)+pDown(1))
	write(20,'(2F15.10)') xi(1),diff
	write(20,*)	pChangeUp,pChangeDown
	xi(1)=xi(1)+dx
	xi(2)=0.0_dp
	pUp(1)= pUp(1) + 10.0_dp*diff
	pDown(1) =  pDown(1) - 10.0_dp*diff
end do

close(19)
close(20)

end subroutine findSep_TBL