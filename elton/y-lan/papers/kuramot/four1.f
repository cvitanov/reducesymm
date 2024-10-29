c This programe is used to give the Fourier componets of sampling points x(1)
c to x(2^n), the frequency increases like ----><----, where the end points 
c are zero frequency and the smallest negative frequency. But this is exactly
c what we want! n has to be a power of 2. You have to divide the result by nn
c to get the inverse FT.


      subroutine four1(nn,isign,dat)
      integer nn,isign,dots,i,j,m,n,istep,mmax 
      double precision dat(2*nn),wtemp,wr,wi,wpr,wpi
      double precision theta,tempr,tempi

      n=nn*2
      j=1
c      write(*,*) 'in',dat
      do 1000,i=1,n-1,2
        if (j.gt.i) then
          wtemp=dat(i)
          dat(i)=dat(j)
          dat(j)=wtemp
          wtemp=dat(i+1)
          dat(i+1)=dat(j+1)
          dat(j+1)=wtemp
        end if
        m=n/2
        do 1100 while(m.ge.2.and.j.gt.m)
          j=j-m
          m=m/2
1100    continue
        j=j+m
1000  continue
c      write(*,*)'middle',dat
c Here begins the danielson-lanczos method.
      mmax=2
      do 2000 while(n.gt.mmax)
        istep=mmax*2
        theta=isign*(6.28318530717959/mmax)
        wtemp=sin(0.5*theta)
        wpr=-2.0*wtemp*wtemp
        wpi=sin(theta)
        wr=1.0
        wi=0.0
        do 2100 m=1,mmax-1,2
          do 2110 i=m,n,istep
            j=i+mmax
            tempr=wr*dat(j)-wi*dat(j+1)
            tempi=wr*dat(j+1)+wi*dat(j)
            dat(j)=dat(i)-tempr
            dat(j+1)=dat(i+1)-tempi
            dat(i)=dat(i)+tempr
            dat(i+1)=dat(i+1)+tempi
2110      continue
          wtemp=wr
          wr=wtemp*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
c          write(*,*)'wr,wi',wr,wi,dat
2100    continue
c        write(*,*) 'end',dat
        mmax=istep
2000  continue
      end











