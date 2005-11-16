module ifc_findroot

interface
	subroutine zbrac(func,x1,x2,succes)
	use nrtype; use nrutil, only : nrerror
	implicit none
	real(sp), intent(inout) :: x1,x2
	logical(lgt), intent(out) :: succes
	interface
		function func(x)
		use nrtype
		implicit none
		real(sp), intent(in) :: x
		real(sp) :: func
		end function func
	end interface
	end subroutine
end interface

end module