PROGRAM FD

USE nrtype
USE interfcs

IMPLICIT NONE

INTEGER(I4B) :: N, d ,i,j,m1,m2,m,ni, dum,indbnd=0, n_iter,max_impr=32, good_impr=10,  updt=1
INTEGER(I4B), ALLOCATABLE :: indx(:)
REAL(DP), ALLOCATABLE ::  LHM(:,:), dumLHM(:,:), vel(:), MVdum(:,:), b(:), b_impr(:),av(:)
REAL(DP), ALLOCATABLE ::  deltaA(:),deltaAn(:),d_deltaAn(:),U(:,:),Usaved(:,:),V(:,:), al(:,:)
REAL(DP) :: lambda, deltatau=0.01_dp, ch, con=0,cons=0, con_old=0, cons_old=0, dtau_update=0.0_dp
REAL(DP) :: t=0, tollerance=1E-6_DP,dd=0.0_dp
REAL(DP) :: iter_tollerance=1E-8_DP, dev=0.0_dp
REAL(DP), ALLOCATABLE :: Ho(:,:),Hinv(:,:), q(:)
INTEGER(I4B) :: just_dec=0, n_no_dec=0, opt_no_dec=4
REAL(DP) :: change_rate=0.02_dp
REAL(DP) :: safe_factor=0.1_dp


N=500
d=16
lambda= 0.8707_dp/(2*PI)
m1=2*d
m2=2*d
ch=N/(12*2*PI)
indx=0
al=0




allocate( LHM(N*d+1,m1+1+m2), dumLHM(N*d+1,m1+1+m2), b(N*d+1), b_impr(N*d+1), vel(d) )

allocate( av(N*d+1), q(N*d+1) )

allocate(deltaA(N*d+1),deltaAn(N*d+1),d_deltaAn(N*d+1),indx(N*d+1), al(N*d+1,m1))

allocate(U(N*d+1,4*d+1),Usaved(N*d+1,4*d+1),V(N*d+1,4*d+1),Ho(4*d+1,4*d+1),Hinv(4*d+1,4*d+1))


av=0
MVdum=0
LHM=0
dumLHM=0
b=0
deltaA=0
deltaAn=0
d_deltaAn=0
U=0
Usaved=0
Ho=0
Hinv=0
V=0
b_impr=0


OPEN(9,FILE='LoopGuess.dat')

	do i=1,N
		Read(9,*)	 av((i-1)*d+1:i*d)
	end do	 

CLOSE(9)

ni=0

call setUV(N,d,ch,U,V)

call setRHS(N,d,m1,m2,av,lambda,U,V,ch,b) ! create the RHS and the correction matrices U and V (without the velocity terms)

!check convergence
con=vector_abs_max(b)

cons=0
do i=1,size(b)
 cons=cons+b(i)**2
end do
cons=sqrt(cons)/N


Print *,ni, "Period", 2*Pi*lambda,"conv1",con,"conv", cons


do while (cons > tollerance)

	ni=ni+1


	b=deltatau*b ! the fictitious timestep was left out


	!add the boundary correction terms (velocities) to U, Usaved and V
	Usaved(:,4*d+1)=0
	do i=1,N-2  
		call derivs(t,av((i-1)*d+1:i*d),vel)
		U((i-1)*d+1:i*d,4*d+1)=-vel
		Usaved((i-1)*d+1:i*d,4*d+1)=-vel
	end do


	V(N*d+1,4*d+1)=1

	! so that now we can form the full Left Hand Side Matrix
	call setLHM(N,d,m1,m2,av,lambda,ch,LHM) 

	updt=1

	if (indbnd == 1) then ! use approximate decomposition from previous step
		n_no_dec=n_no_dec+1
		if ( mod(n_no_dec,30) == 0 ) then 
			dtau_update=(1+(n_no_dec-opt_no_dec)*change_rate)*deltatau
			Print *, "updated deltatau because of good convergence to", deltatau
		end if
		!call banbks(dumLHM,m1,m2,al,indx,b_impr) ! backsubstitute with the last LU decomposition
		call fastWoodbury(dumLHM,Usaved,V,Hinv,N,d,m1,m2,deltaA,al,indx,deltaAn) ! Apply woodbury to find deltaAn (here n=0)
		call banmul(LHM,m1,m2,deltaAn,b_impr)	! to check approximation
		b_impr = b - (deltaA+MatMul( U , MatMul(Transpose(V),deltaA) )  )  ! to check approximation and also form the next rhs
		dev=vector_abs_max(b_impr)	! measure of approximation
		!Print *, "impr", dev
		n_iter=0
		do while( (dev > iter_tollerance) .AND. (n_iter<max_impr) ) ! iterate until good approximation or too many iterations
			n_iter=n_iter+1	! count iterations
			call banbks(dumLHM,m1,m2,al,indx,b_impr)	! backsubstitute with the last LU decomposition
			call instantWoodbury(Usaved,V,Hinv,b_impr,d_deltaAn) ! Find correction to deltaAn due to Woodbury
			deltaAn=deltaAn+d_deltaAn ! and apply it
			call banmul(LHM,m1,m2,deltaAn,b_impr) ! to check approximation
			b_impr=b - ( b_impr+MatMul(U,MatMul(Transpose(V),deltaAn)) )	! to check approximation and also form the next rhs
			dev=vector_abs_max(b_impr)	! measure of approximation
		!	Print *, "iter impr", dev
		end do
		if ( n_iter>=good_impr ) then !if number of iterations exceeds tolerable value without (or even with) good convergence
			if ( dev < iter_tollerance ) then
				Print *, " convergence after ", n_iter, "iterations."
				updt=1 ! keep result
			else
				Print *, " NO convergence after ", n_iter, "iterations." 
				updt=0	! Don't keep result
			end if
			indbnd=0 ! go for a fresh LU decomposition
		end if
	else
		dtau_update=(1+(n_no_dec-opt_no_dec)*change_rate)*deltatau
		if ( ( dtau_update< 1.0_dp )  .and. ( dtau_update > 0.0_dp ) ) then
			deltatau=dtau_update  ! adaptive timestep
		end if
		Print *, "Decomposed after", n_no_dec, "approximate inversions."
		Print *, "Delta tau=", deltatau 
		n_no_dec=0 
		dumLHM=LHM
		call bandec(dumLHM,m1,m2,al,indx,dd) ! Performs the banded LU Decomposition
		!b_impr=b !save it for comparison purposes	
		!Solve the auxiliary problem LHM.y=b
		call banbks(dumLHM,m1,m2,al,indx,b)
		Usaved=U 
		Print *, "before Woodbury"
    	call Woodbury(dumLHM,Usaved,V,N,d,m1,m2,b,al,indx,Ho,Hinv,deltaA) ! applies Woodbury formula
		Print *, "after Woodbury"
		!call banmul(LHM,m1,m2,deltaA,bdum) 
		!Print *, "imprLU", vector_abs_max(bdum+MatMul(U,MatMul(Transpose(V),deltaA))-b)
		indbnd=1
		just_dec=1
	end if

	if (updt==1) then
		av(1:N*d)=av(1:N*d)+deltaA(1:N*d)	! update av
		lambda=lambda+deltaA(N*d+1)			! update lambda	
	else
		Print *, "NO UPDATE"	 			
	end if

	call setUV(N,d,ch,U,V)

	call setRHS(N,d,m1,m2,av,lambda,U,V,ch,b) ! create the RHS and the correction matrices U and V (without the velocity terms)

	con_old=con
	cons_old=cons

	!check convergence
	con=vector_abs_max(b)

	cons=0
	do i=1,size(b)
	 cons=cons+b(i)**2
	end do

	cons=sqrt(cons)/N

	Print *,ni, "Period", 2*Pi*lambda,"conv1",con,"conv", cons

	if (  1.0_dp*cons_old<cons ) then	!revert to old solution
		Print *, "!!!!!!!!!!WARNING!!!!!!!!!!!!"
		OPEN (9,FILE='beforebadcycle.dat') !
			do i=1,N
				 write(9,"(16F12.8)") av((i-1)*d+1:i*d)
			end do
		CLOSE (9)
		!export the derivative at each point
		call banmul(FDM(N,d,m1,m2,ch),m1,m2,av,q)
		call setUV(N,d,ch,U,V)
		q=q + MatMul(U,MatMul(Transpose(V),av))
		OPEN (9,FILE='beforebadcycleD.dat') !
			do i=1,N
				 write(9,"(16F12.8)") q((i-1)*d+1:i*d)
			end do
		CLOSE (9)
		deltatau=deltatau*safe_factor
		av(1:N*d)=av(1:N*d)-deltaA(1:N*d)	!unupdate av
		lambda=lambda-deltaA(N*d+1)			!unupdate lambda	
		if ( (indbnd==1) .and. (just_dec==0) ) then
			call setUV(N,d,ch,U,V)
			call setRHS(N,d,m1,m2,av,lambda,U,V,ch,b) ! create the RHS and the correction matrices U and V (without the velocity terms)
			con=vector_abs_max(b)
			cons=0
			do i=1,size(b)
				cons=cons+b(i)**2
			end do
			cons=sqrt(cons)/N
			indbnd=0
		else
			Print *,"failed with convergence 1", con ,"conv all", cons 
			OPEN (9,FILE='badcycle.dat') !
				do i=1,N
					 write(9,"(16F12.8)") av((i-1)*d+1:i*d)
				end do
			CLOSE (9)
			!export the derivative at each point
			call banmul(FDM(N,d,m1,m2,ch),m1,m2,av,q)
			call setUV(N,d,ch,U,V)
			q=q + MatMul(U,MatMul(Transpose(V),av))
			OPEN (9,FILE='badcycleDerivative.dat') !
				do i=1,N
					 write(9,"(16F12.8)") q((i-1)*d+1:i*d)
				end do
			CLOSE (9)			 	 
			call exit
		end if	 			
	end if
	
	just_dec=0

end do

OPEN (9,FILE='cycle.dat')

	do i=1,N
		 write(9,"(16F12.8)") av((i-1)*d+1:i*d)
	end do

CLOSE (9)

END PROGRAM