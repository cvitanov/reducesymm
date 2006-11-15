subroutine SetANdiag_KS_diff(a,Andiag)

use nrtype
use ifc_int_ks
use nrutil, only:assert_eq

implicit none

complex(dpc), dimension(:), intent(in) :: a
real(dp), dimension(:,:), intent(out) :: Andiag
!
!

complex(dp), dimension(size(a)+1) :: adum,N_adum_p,N_adum_m
integer(i4b):: k,j,d,i 
real(dp) :: tmp,da=0.0000001_dp


d=assert_eq(size(Andiag,1),size(Andiag,2),2*size(a),'SetANdiag')

adum=(0.0_dp,0.0_dp)

Andiag=0.0_dp

do i=1,d/2	
	adum(2:d/2)=a
	tmp=real(a(i))+da
	call donothing(tmp)
	da=tmp-real(a(i))
	adum(i+1)=real(a(i))+da + ii*aimag(a(i)) 
	call SetNlin_KS(adum,N_adum_p) 
	adum(i+1)=real(a(i))-da + ii*aimag(a(i))
	call SetNlin_KS(adum,N_adum_m) 
	do k=1,d/2
		Andiag(k,i)=(real(N_adum_p(k+1))-real(N_adum_m(k+1)))/(2.0_dp*da)
	end do
	do k=1,d/2	
		Andiag(d/2+k,i)=(aimag(N_adum_p(k+1))-aimag(N_adum_m(k+1)))/(2.0_dp*da)
	end do
end do

do i=1,d/2	
	adum(2:d/2)=a
	tmp=aimag(a(i))+da
	call donothing(tmp)
	da=tmp-aimag(a(i))
	adum(i+1)=real(a(i))+ii*(aimag(a(i))+da)
	call SetNlin_KS(adum,N_adum_p) 
	adum(i+1)=real(a(i))+ii*(aimag(a(i))-da) 
	call SetNlin_KS(adum,N_adum_m) 
	do k=1,d/2
		Andiag(k,d/2+i)=(real(N_adum_p(k+1))-real(N_adum_m(k+1)))/(2.0_dp*da)
	end do
	do k=1,d/2	
		Andiag(d/2+k,d/2+i)=(aimag(N_adum_p(k+1))-aimag(N_adum_m(k+1)))/(2.0_dp*da)
	end do
end do

end subroutine
