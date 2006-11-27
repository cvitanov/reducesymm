program integrManif

use nrtype
use ifc_integr
use f95_lapack, only: LA_GEEV
use la_precision, only: wp => dp
use ifc_util
use ks

implicit none

include "fftw3.f"

real(dp), dimension(:), allocatable :: v, vdum, bc, bc0,bc3,fvec,bcdum
complex(dpc), dimension(:),allocatable :: a,adum,w
complex(dpc), dimension(:), allocatable :: ai,af
complex(dpc), dimension(:,:), allocatable :: jacdum, VR, VR2
integer(i8b) :: invplan, plan ! needed by fftw3
integer(i4b) :: k,i,j,q,sdim,Nrep, Np 
real(dp) :: T,h,ti,tf,dst,dst0
character*64 :: wd
integer(i4b) :: nargs
logical :: logicdum
integer(i4b) :: stbl_dir,unst_dir

nargs=iargc()

if (nargs .ne. 1) then
	print *,"Program must be called with exactly one argument indicating the data directory."
	call exit
end if

call getarg(1,wd)

open(21,file=trim(wd)//'/parameters.dat')
	read(21,*)
	read(21,*) d
	read(21,*)
	read(21,*) L 
	read(21,*)
	read(21,*) Nsteps 
	read(21,*)
	read(21,*) T
	read(21,*)
	read(21,*) Nrep 
	read(21,*) 
	read(21,*) Nplt 
	read(21,*)
	read(21,*) Np 
	read(21,*)
	read(21,*) dst0 
	read(21,*)
	read(21,*) dst  
	read(21,*)
	read(21,*) M
	read(21,*)
	read(21,*) R
close(21)

212 Format(<d>F30.18)
211 Format(<d/2+1>F30.18)
230 Format(<q>F30.18)

allocate(v(d),vdum(d),bc(d),bc0(d),bc3(d),fvec(d),bcdum(d))
allocate(a(d/2+1),adum(d/2+1),ai(d/2+1),af(d/2+1))
allocate(lin(d/2+1),f0(d/2+1),f1(d/2+1),f2(d/2+1),f3(d/2+1),e(d/2+1),e2(d/2+1))
allocate(w(d),vR(d,d),vR2(d,d),jac(d,d),jacdum(d,d))

a=0.0_dp

open(19,file=trim(wd)//'Uequil-2w.dat')
 
	read(19,*) v(1:d)
 
close(19)

tf=real(Nrep,dp)*T
ti=0.0_dp
h=abs(tf-ti)/Nsteps

call dfftw_plan_dft_r2c_1d(plan,d,v,a,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
a=a/size(v)

bc0(1:d/2)=real(a(2:size(a)))
bc0(d/2+1:d)= aimag(a(2:size(a)))

call dfftw_plan_dft_r2c_1d(plan,d,v,a,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
a=a/size(v)

bc3(1:d/2)=real(a(2:size(a)))
bc3(d/2+1:d)= aimag(a(2:size(a)))


call ksFJ(bc0,fvec,jac)

jacdum=jac

call la_geev(jacdum,w,VR=VR)

print *,"w0",w(119),w(120),w(128),w(125)

VR2=VR

call sort_pick_2Re(w,VR)

print *,"w",w(120:128)

print *,"dot",dot_product(VR2(:,120),VR(:,127))
print *,"dot",dot_product(VR2(:,125),VR(:,123))

stbl_dir=122
unst_dir=127

open(15,file=trim(wd)//'eigvalues-2w.dat')
	do i=1,d
		write(15,*) w(i)
	end do
close(15)

! open(29,file='eigvector-2w.dat')
! do i=1,
! 	write(29,*) vR(:,127)
! enddo
! close(29)

call SetLin_KS(Lin)

call etdrk4DiagPrefactors(Lin,h,R,M,f0,f1,f2,f3,e,e2)

open(17,file=trim(wd)//'aP.dat')
open(18,file=trim(wd)//'a123.dat')

do i=1,Np
	print *,"ic",i
	bc(:)=bc0+real(dst0+(i-1)*dst,dp)*(real(vR(:,unst_dir))+aimag(vR(:,unst_dir)))
	
	ai(2:size(a))=bc(1:d/2)+ii*bc(d/2+1:d)

	call etdrk4DiagDriverS(ti,ai,Nsteps,tf,af,f0,f1,f2,f3,e,e2,Nplt,SetNlin_KS)
	do k=1,size(aSt,1)
		q=8
		write(18,230) real(aSt(k,2)),aimag(aSt(k,2)),real(aSt(k,3)),aimag(aSt(k,3)),real(aSt(k,4)),aimag(aSt(k,4)),real(aSt(k,5)),aimag(aSt(k,5))
		bcdum(1:d/2)=real(aSt(k,2:size(a)))
		bcdum(d/2+1:d)= aimag(aSt(k,2:size(a)))
		q=3
		write(17,230)  dot_product(real(VR(:,unst_dir)),bcdum),dot_product(aimag(VR(:,unst_dir)),bcdum),dot_product(real(VR(:,stbl_dir)),bcdum)
	end do
end do


close(17)
close(18)


end program