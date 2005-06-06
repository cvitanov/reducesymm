PROGRAM integr

USE nrtype
USE ifc_integr, ONLY: rk4Pdriver, rk4driver, rk2Jdriver
USE ifc_util, ONLY: UnitMatrix
USE parameters

IMPLICIT NONE

REAL(DP), DIMENSION(:,:), ALLOCATABLE :: y, yclose, y_poinc,  J
REAL(DP) :: xi=0.0_dp, xf=0.1_dp, xf_poinc, jxi,jxf
INTEGER(I4B) :: nsteps=1000, nsteps_poinc=1000, jnsteps=10000, i, sect, n_inters=0, total_inters=2
REAL(DP), DIMENSION(:), ALLOCATABLE :: jyi
INTERFACE
	SUBROUTINE roesslerField(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE roesslerField
END INTERFACE

INTERFACE
	SUBROUTINE roesslerFieldP(x,y,dydx,kappa)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x,kappa
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(:), INTENT(OUT) :: dydx	
	END SUBROUTINE roesslerFieldP
END INTERFACE

interface
	function roesslerVar(x,y)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP), DIMENSION(:), INTENT(IN) :: y
		REAL(DP), DIMENSION(size(y),size(y)) :: roesslerVar
	end function 
end interface



INTEGER(I4B) :: idum, p=0

allocate(y(nsteps+1,d+1), yclose(jnsteps+1,d))
allocate(y_poinc(nsteps_poinc+1,d+1))
allocate(jyi(d+1))
allocate(J(d,d))

y=0.0_dp

open(9,file='initial.dat')

read(9,*)  y(1,1:d)

close(9)

y(1,d+1)=xi

sect=1

open (10, file='section.dat')
	read(10,*)	xf_poinc
close(10)

!OPEN(10, file='roessler.dat')
!OPEN(11, file='roessler-poinc.dat')

do while( (n_inters<total_inters) )
	call rk4driver(xi,y(1,:),xi+xf,nsteps,y,roesslerField)
	if  ( ( y(1,sect) > xf_poinc  ) .AND. ( y(nsteps,sect) < xf_poinc ) ) THEN
		p=1
		call rk4Pdriver(y(1,sect),y(1,:),xf_poinc,nsteps_poinc,y_poinc,roesslerFieldP,p,sect)
! 		DO i=1, nsteps_poinc
! 			WRITE(10,'(4F13.7)') y_poinc(i,:) 
! 		END DO
		n_inters=n_inters+1
		p=0
		Write (11,'(4F13.7)') y_poinc(size(y_poinc,1),:)
		y=0.0
		y(1,:)=y_poinc(size(y_poinc,1),:)
		y(1,sect)=xf_poinc
		if (n_inters==1) then
			jxi = y_poinc(size(y_poinc,1),d+1)
			jyi = y_poinc(size(y_poinc,1),:)
		else if (n_inters==2) then
			jxf = y_poinc(size(y_poinc,1),d+1)
		end if
		xi=xi+xf
		y_poinc=0
!		Read *,idum
	else
! 		DO i=1, nsteps
! 			WRITE(10,'(4F13.4)') y(i,:) 
! 		END DO
		y(1,:)=y(size(y,1),:)
		y(2:size(y,1),:)=0
		xi=xi+xf
		Print *,xi
	end if
end do

J=UnitMatrix(d)


deallocate(y)
allocate(y(jnsteps+1,d+1))

p=0

call rk2Jdriver(jxi,jyi,jxf,jnsteps,y,J,J,roesslerVar,roesslerField)

open(7,file="cycle.dat")
do i=1,size(y,1)
	write(7,format_label) y(i,1:d)
end do
close(7)

open(8,file="pointPoinc.dat")
	write(8,format_label) y(1,1:d)
close(8)

open(9,file='J.dat')

do i=1,d
	write(9,format_label) J(i,:)
end do

close(9)

!CLOSE(11)

!CLOSE(10) 


END PROGRAM