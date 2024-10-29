c The routine 'dftint' is used to calculate the Fourier integral of a function
c represented by dat1(.), giving the real and imaginary parts of the resulting
c wavenumber components. The routine 'dftcor' is used to calculate the 
c corrections for computing 
c Fourier integrals of periodic functions utilizing FFT.

      subroutine dftint(n,dat1,dat2,datk1,datk2)
      integer n,i,j,isign,m1,m2
      double precision dat1(n),dat2(n),datk1(2*n)
      double precision datk2(2*n),endpts2(8) 
      double precision a0,b0,pi,omg,delta,endpts1(8)
      double precision t1,t2,t3,t4,crel,cimg,corfac

      pi=3.1415926535897932D0
      a0=0
      b0=2*pi
      delta=(b0-a0)/n 
      do 1000 i=1,n,1
        j=2*(i-1)
        datk1(j+1)=dat1(i)
        datk1(j+2)=dat2(i)
1000  continue
      isign=1
      call four1(n,isign,datk1)
      m1=2*n+2
      m2=m1+1
      datk2(1)=datk1(2)
      datk2(2)=0 
      datk1(2)=0
      do 2000 j=3,n+1,2
        t1=0.5*(datk1(j)+datk1(m1-j))
        t2=0.5*(datk1(j+1)-datk1(m2-j))
        t3=0.5*(datk1(j+1)+datk1(m2-j))
        t4=-0.5*(datk1(j)-datk1(m1-j))
        datk1(j)=t1
        datk1(j+1)=t2
        datk2(j)=t3
        datk2(j+1)=t4
        datk1(m1-j)=t1
        datk1(m2-j)=-t2
        datk2(m1-j)=t3
        datk2(m2-j)=-t4     
2000  continue
      do 3000 i=1,4,1
        endpts1(i)=dat1(i)
        endpts2(i)=dat2(i)
3000  continue
      do 4000 i=1,3,1 
        endpts1(i+4)=dat1(n-3+i)
        endpts2(i+4)=dat2(n-3+i)
4000  continue
      endpts1(8)=dat1(1)
      endpts2(8)=dat2(1)
      t1=2*pi
      do 5000 i=0,n/2-1,1
        omg=i
        m1=2*i
        call dftcor(omg,delta,a0,b0,endpts1,crel,cimg,corfac) 
        dat1(m1+1)=delta*(corfac*datk1(m1+1)+crel)/t1
        dat1(m1+2)=-delta*(corfac*datk1(m1+2)+cimg)/t1
        call dftcor(omg,delta,a0,b0,endpts2,crel,cimg,corfac) 
        dat2(m1+1)=delta*(corfac*datk2(m1+1)+crel)/t1
        dat2(m1+2)=-delta*(corfac*datk2(m1+2)+cimg)/t1
5000  continue
      end

      subroutine dftcor(omg,delta,a,b,endpts,crel,cimg,corfac)
      double precision omg,a,b,crel,cimg,corfac,delta
      double precision endpts(8),sth4i,tth4i,tht
      double precision t1,t2,t3,t4,t5,t6,tmt2,spt2
      double precision a0i,a0r,a1i,a1r,a2i,a2r,a3i,a3r
      double precision cl,sl,cr,sr,arg,athr,athi 

      tht=omg*delta
      if (a.ge.b.or.(tht.lt.0.or.tht.gt.3.1416)) then
        write(*,*) 'Bad arguments to dftcor'
        stop
      end if
      if (abs(tht).lt.5.0e-2) then
        t2=tht*tht
        t4=t2*t2
        t6=t4*t2 
        corfac=1.0-(dble(11)/720.0)*t4+(dble(23)/15120.0)*t6
        a0r=(-dble(2)/3.0)+t2/45.0+(dble(103)/15120.0)*t4
     \-(dble(169)/226800.0)*t6
        a1r=(dble(7)/24)-(dble(7)/180.0)*t2+(dble(5)/3456.0)*t4
     \-(dble(7)/259200.0)*t6
        a2r=(-dble(1)/6.0)+t2/45.0-(dble(5)/6048.0)*t4+t6/64800.0
        a3r=(dble(1)/24.0)-t2/180.0+(dble(5)/24192.0)*t4-t6/259200.0
        a0i=tht*((dble(2)/45.0)+(dble(2)/105.0)*t2-(dble(8)/2835.0)*t4
     \+(dble(86)/467775.0)*t6)
        a1i=tht*((dble(7)/72.0)-t2/168.0+(dble(11)/72576.0)*t4
     \-(dble(13)/5987520.0)*t6)
        a2i=tht*(-(dble(7)/90.0)+t2/210.0-(dble(11)/90720.0)*t4
     \+(dble(13)/7484400.0)*t6)
        a3i=tht*((dble(7)/360)-t2/840.0+(dble(11)/362880.0)*t4
     \-(dble(13)/29937600.0)*t6) 
        athr=dble(1)/3.+t2/45.-8.*t4/945.+11*t6/14175.
        athi=-a0i
      else
        t1=cos(tht)
        t3=sin(tht)
        t5=t1**2-t3**2
        t6=2*t1*t3
        t2=tht*tht
        t4=t2*t2
        tmt2=3-t2
        spt2=6+t2
        sth4i=1.0/(6*t4)
        tth4i=2*sth4i
        corfac=tth4i*spt2*(3-4*t1+t5)
        a0r=sth4i*(-42+5*t2+spt2*(8*t1-t5))
        a0i=sth4i*(tht*(-12+6*t2)+spt2*t6)
        a1r=sth4i*(14*tmt2-7*spt2*t1)
        a1i=sth4i*(30*tht-5*spt2*t3)
        a2r=tth4i*(-4*tmt2+2*spt2*t1)
        a2i=tth4i*(-12*tht+2*spt2*t3)
        a3r=sth4i*(2*tmt2-spt2*t1)
        a3i=sth4i*(6*tht-spt2*t3) 
        athr=sth4i*(-6+11*t2+spt2*t5)
        athi=-a0i
      end if 
      cl=a0r*endpts(1)+a1r*endpts(2)+a2r*endpts(3)+a3r*endpts(4)
      sl=a0i*endpts(1)+a1i*endpts(2)+a2i*endpts(3)+a3i*endpts(4)
      cr=athr*endpts(8)+a1r*endpts(7)+a2r*endpts(6)+a3r*endpts(5)
      sr=athi*endpts(8)-a1i*endpts(7)-a2i*endpts(6)-a3i*endpts(5)
      arg=omg*(b-a)
      t1=cos(arg)
      t2=sin(arg)
      crel=cl+t1*cr-t2*sr
      cimg=sl+t2*cr+t1*sr      
      end 
