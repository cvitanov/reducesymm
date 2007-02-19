	SUBROUTINE mnewt_elim(ntrial,elim,x,tolx,tolf,usrfun)
	
	USE nrtype
	use ifc_newt
	use la_precision, only: wp => dp
	use f95_lapack, only: LA_GESV

	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: ntrial
	integer(i4b), DIMENSION(:), INTENT(IN) :: elim ! variable to eliminate (not set up to use it yet)
	REAL(DP), INTENT(IN) :: tolx,tolf
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: x

	INTERFACE
		SUBROUTINE usrfun(x,elim,fvec,fjac)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		integer(i4b), DIMENSION(:), INTENT(IN) :: elim 
		REAL(DP), DIMENSION(:), INTENT(OUT) :: fvec
		REAL(DP), DIMENSION(:,:), INTENT(OUT) :: fjac
		END SUBROUTINE usrfun
	END INTERFACE

	INTEGER(I4B) :: i
	REAL(DP), DIMENSION(size(x)-1) :: fvec,p ! eliminate one variable
	REAL(DP), DIMENSION(size(x)-1,size(x)-1) :: fjac ! eliminate one variable
	real(dp) :: sumx,sumf
	
	newton_condition_met=0

	do  i=1,ntrial
		call usrfun(x,elim,fvec,fjac)
		sumf=sum(abs(fvec))
		if ( sumf <= tolf) then
			print *,"Condition for fvec met, with fvec=",sumf
			newton_condition_met=1
			RETURN
		endif
		p=-fvec
		call la_gesv(fjac,p)
		x(2:size(x))=x(2:size(x))+p
		sumx=sum(abs(p))
		print *,"Delta x",sumx
		if (sumx <= tolx) then 
			print *,"Condition for tolx met, with tolx=",sumx
			newton_condition_met=2
			RETURN
		end if
	end do

	END SUBROUTINE mnewt_elim
