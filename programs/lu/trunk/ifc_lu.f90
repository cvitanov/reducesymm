module ifc_lu

INTERFACE
	SUBROUTINE lubksb(a,indx,b)
	USE nrtype
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
	END SUBROUTINE lubksb
END INTERFACE

INTERFACE
	SUBROUTINE ludcmp(a,indx,d)
	USE nrtype
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
	REAL(DP), INTENT(OUT) :: d
	END SUBROUTINE ludcmp
END INTERFACE


end module