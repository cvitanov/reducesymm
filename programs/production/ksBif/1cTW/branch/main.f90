Program ksEquil

use nrtype
use ks
use ifc_newt
use ifc_util
use f95_lapack, only: LA_GEEV
USE LA_PRECISION, ONLY: WP => DP

implicit none

include "fftw3.f"

real(dp), dimension(:),allocatable :: v,vx,vxx,vxxx
complex(dpc), dimension(:),allocatable :: a,adum
complex(dp), dimension(:), allocatable :: w
real(dp), dimension(:),allocatable :: bc,fvec
real(dp), dimension(:,:),allocatable :: fjac
complex(dp), dimension(:,:), allocatable :: fjacdum,vR
real(dp), dimension(:),allocatable :: ar,ai
integer(i8b) :: invplan, plan ! needed by fftw3
integer(i4b) :: k,i,j,sdim=8, Nf
real(dp) :: dL, kappa
character*64 :: wd
integer(i4b) :: nargs
logical :: logicdum

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
	read(21,*) dL
	read(21,*)
	read(21,*) tolbc
	read(21,*)
	read(21,*) tolf
	read(21,*)
	read(21,*) Ntrial
	read(21,*)
	read(21,*) Nsteps
	read(21,*)
	read(21,*) Nplt
	read(21,*)
	read(21,*) M
	read(21,*)
	read(21,*) R
	read(21,*)
	read(21,*) Nf
close(21)

220 Format(F21.16)
221 Format(<d>F21.16)
222 Format(<d/2+1>F21.16)
! 

allocate(v(d),vx(d),vxx(d),vxxx(d),a(d/2+1),adum(d/2+1),bc(d),ar(d/2+1),ai(d/2+1),w(d))
allocate( fvec(d),fjac(d,d),fjacdum(d,d), vR(d,d))
! allocate(v(d),vdum(d),bc(d))
!allocate(a(d/2+1),adum(d/2+1),ai(d/2+1),af(d/2+1))
allocate(lin(d/2+1),f0(d/2+1),f1(d/2+1),f2(d/2+1),f3(d/2+1),e(d/2+1),e2(d/2+1))
allocate(f0dum(d/2+1),f1dum(d/2+1),f2dum(d/2+1),f3dum(d/2+1),edum(d/2+1),e2dum(d/2+1))

open(19,file=trim(wd)//'/ic.dat')
 
	read(19,*) v
 
close(19)

open(20,file=trim(wd)//'/periodsGuess.dat')

        read(20,*) Tw
        read(20,*) kappa

close(20)

call dfftw_plan_dft_r2c_1d(plan,d,v,a,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)
a=a/size(v)

bc(1:d/2)=real(a(2:size(a)))
bc(d/2+1:d)= aimag(a(2:size(a)))

open(35,file=trim(wd)//'/Jdiag.dat')
open(24,file=trim(wd)//'/equilU.dat')
open(26,file=trim(wd)//'/c_steady.dat')
open(27,file=trim(wd)//'/periods.dat')

print *, "T",Tw,"Nsteps",Nsteps,"Nplt",Nplt

do i=1,Nf

	print *, "Point #",i, "L=",L
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call mnewtTW(ntrial,bc,tolbc,tolf,kappa,ksFJ_req)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	print *,"T",Tw,"kappa",kappa
	

	write (27,220) Tw
	write (27,220) kappa


! 	call ksFJ(bc,fvec,fjac)
! 
! 	fjacdum=fjac
! 	w=0.0_dp
! 	vR=0.0_dp
! 	call la_geev(fjacdum,w,vR=vR) ! ,select=SelectSmallEig_c,sdim=sdim)
! 	call sort_pick(w)
! 
! 	write(35,"(4F28.16)") real(w(d)),real(w(d-1)),real(w(d-2)),real(w(d-3))

	a=(0,0)

	a(2:size(a))=bc(1:d/2)+ii*bc(d/2+1:d)

	call dfftw_plan_dft_c2r_1d(invplan,d,a,v,FFTW_ESTIMATE)
	call dfftw_execute(invplan)
	call dfftw_destroy_plan(invplan)

	write(24,221) v
	
	call FourierDif(v,vx,L,1)
	call FourierDif(v,vxx,L,2)
	call FourierDif(v,vxxx,L,3)
	
	write(26,"(2F21.16)") L, sum(v**2-vx-vxxx)/d
	print *,"E",sum(v**2-vx-vxxx)/d
	
	L=L+dL
	
end do

close(35)
close(26)
close(24)
close(27)

end program
