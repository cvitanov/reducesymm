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
integer(i4b) :: k,j,idum,q,sdim,Nrep 
integer(i8b) :: i,No, Nf, Ndiv 
real(dp) :: T,h,ti,tf,dst,dst0
character*64 :: wd
integer(i4b) :: nargs
logical :: logicdum
integer(i4b) :: stbl_dir,unst_dir, direct

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
	read(21,*) No
	read(21,*)
	read(21,*) Nf
	read(21,*)
	read(21,*) Ndiv 
	read(21,*)
	read(21,*) dst0 
	read(21,*)
	read(21,*) dst  
	read(21,*)
	read(21,*) Mi
	read(21,*)
	read(21,*) R
	read(21,*)
	read(21,*) unst_dir
close(21)

212 Format(<d>F30.18)
211 Format(<d/2+1>F30.18)
230 Format(<q>F30.18)

allocate(v(d),vdum(d),bc(d),bc0(d),bc3(d),fvec(d),bcdum(d))
allocate(a(d/2+1),adum(d/2+1),ai(d/2+1),af(d/2+1))
allocate(lin(d/2+1),f0(d/2+1),f1(d/2+1),f2(d/2+1),f3(d/2+1),e(d/2+1),e2(d/2+1))
allocate(w(d),vR(d,d),vR2(d,d),jac(d,d),jacdum(d,d))

a=0.0_dp

open(19,file=trim(wd)//'equilU.dat')
 
	read(19,*) v(1:d)
 
close(19)

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

tf=real(Nrep,dp)*TWOPI_D/aimag(w(unst_dir))
ti=0.0_dp
h=abs(tf-ti)/Nsteps
print *,"timestep", h

call SetLin_KS(Lin)

call etdrk4DiagPrefactors(Lin,h,R,Mi,f0,f1,f2,f3,e,e2)

print *,"distance",dst0,dst

do idum=1,100
	open(17,file=trim(wd)//'aP.dat')
	open(18,file=trim(wd)//'a123.dat')
	open(19,file=trim(wd)//'v.dat')
	do i=No,Nf
		print *,"No=",No,"Nf=",Nf,"Ndiv=",Ndiv
		print *,"ic",i,real(w(unst_dir))*(TWOPI_D/abs(aimag(w(unst_dir))) )*(i-1)/Ndiv
		bc(:)=bc0+dst0*Exp(real(w(unst_dir))*(TWOPI_D/abs(aimag(w(unst_dir))) )*(i-1)/Ndiv)*real(vR(:,unst_dir))
	!	bc(:)=bc0+dst0*Exp(3.66343*(i-1)/Np)*real(vR(:,unst_dir))  !+aimag(vR(:,unst_dir)))
	!       bc(:)=bc0+dst0*( cos(2*PI*(i-1)/Np)*real(vR(:,unst_dir)) + sin(2*PI*(i-1)/Np)*aimag(vR(:,unst_dir)) )
		
		ai(2:size(a))=bc(1:d/2)+ii*bc(d/2+1:d)
	
		call etdrk4DiagDriverS(ti,ai,Nsteps,tf,af,f0,f1,f2,f3,e,e2,Nplt,SetNlin_KS)
		do k=1,size(aSt,1)
			q=8
			write(18,230) real(aSt(k,2)),aimag(aSt(k,2)),real(aSt(k,3)),aimag(aSt(k,3)),real(aSt(k,4)),aimag(aSt(k,4)),real(aSt(k,5)),aimag(aSt(k,5))
			bcdum(1:d/2)=real(aSt(k,2:size(a)))
			bcdum(d/2+1:d)= aimag(aSt(k,2:size(a)))
			q=3
			write(17,230)  dot_product(real(VR(:,unst_dir)),bcdum),dot_product(aimag(VR(:,unst_dir)),bcdum),dot_product(real(VR(:,stbl_dir)),bcdum)
			adum=aSt(k,:)
			call dfftw_plan_dft_c2r_1d(invplan,d,adum,v,FFTW_ESTIMATE)
			call dfftw_execute(invplan)
			call dfftw_destroy_plan(invplan)
			write(19,"(<d>F20.16)") v
		end do
	end do
	close(17)
	close(18)
	close(19)
	print *,"Direction (0,1)?"
	read(*,*) direct
	Ndiv=2*Ndiv
	No=No+direct
	No=2*No-1
	Nf=No+2
end do



end program
