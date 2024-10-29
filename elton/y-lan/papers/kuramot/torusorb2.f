c This progam employs the first order partial implicit method to evolve
c the orbit in the phase space.  The requirement is that
c d,even, trunc,2^n  

      program torusorb2
      implicit none
      integer i,j,k,d,trunc,tp,dtr,flag,count1,m,nn,count2,tr
      integer kdtr,count3
      parameter(d=4,trunc=32,nn=2,tr=trunc**nn,dtr=d*tr,m=1)
      parameter(kdtr=d*(trunc/2-1)**2)
      double precision tiny,hmin,eps,dh,hdid,hnext
      double precision nu,dydx(d),index(kdtr+nn)
      double precision a(d),x(dtr),xx(dtr),time,time1
      integer m1,m2,m3,resterm,count,ct2,pn,cycle
      double precision p,pi,jb(d,d,tr)
      double precision dv(d,d),limit,invtr,y(dtr)
      double precision gtime,t4,t44,t0,x0(dtr/m),dat2(2*trunc/m)
      double precision x01(trunc/m),x02(trunc/m),dat1(2*trunc/m)
      double precision eps1,eps2,eps3,omg1,omg2
      double precision matr(kdtr+nn,kdtr+nn)

      data tiny,hmin,eps/1.0e-30,1.0e-20,1e-4/
c  eps is used to control the tolerence error.
  
      pi=3.1415926535897932D0
      gtime=2*pi
c      dh=1e-3
      nu=-0.015
      cycle=13
      resterm=2
      eps1=0.07
      eps2=0.1
      eps3=0.004
      limit=.5D-4
      invtr=dble(1.)/trunc
      t0=2*pi/trunc
      omg1=2*pi*((sqrt(5.D0)-1.)/2.)
      omg2=2*pi*((sqrt(3.D0)-1.)/2.)
     

c*************** This part is used to initialize the orbit***************
c***Fourier transform the evolution orbit then erase the high-f part **** 
c      open(2,file='antorb4b.dat',status='old')
c      open(3,file='torus22b.dat',status='old')
      open(4,file='tor4d2.dat',status='unknown') 
      do 500 i=0,trunc-1,1
        m1=d*trunc*i
        t4=0.000+i*t0+0.03*sin(i*t0)
        t44=omg1+0.04*cos(i*t0)
        do 510 j=0,trunc-1,1
          m2=m1+d*j
          x(m2+1)=t4
          x(m2+2)=t44
          x(m2+3)=0.000+j*t0+0.0*sin(j*t0)
          x(m2+4)=omg2+0.0*cos(j*t0) 
510     continue
500   continue
c      do 300 i=1,(dtr/m+1)*35,1
c        read(3,'(F12.6)') t1
c300   continue
c      read(3,'(F12.6)') t1
c      do 400 i=1,dtr/m,1
c        read(3,'(F12.6)') x0(i)
c400   continue
c      do 600 i=1,trunc/m,1
c        x01(i)=x0(2*i-1)-(i-1)*t0*m
c        x02(i)=x0(2*i)
c600   continue
CCCCCNote the dimension of x and xx, not general!!!!!
c      call interpfour(trunc/m,m,x01,x02,x,xx,dat1,dat2)
c      do 700 i=1,trunc,1
c        x(2*i-1)=x(2*i-1)+(i-1)*t0
c700   continue
c      close(3,status='keep')
c      write(*,*) 'before',x0(1),x0(2),x0(3),x0(4)
c      write(*,*) 'after',x(1),x(2),x(5),x(6)
c      stop
      go to 20
CCCCCFinish matching time series and start to initialize the dataCCCCC



c      write(*,*) 'before',(x(i),i=1,9*d+1,d)
c      call inismooth(d,trunc,x)
c      do 3000 i=1,(d+1)/2,1
c          m1=2*(i-1)
c        do 3100 j=1,trunc,1
c          m2=2*(j-1)
c          m3=m1+(j-1)*d
c          temp1(m2+1)=x(m3+1)
c          if (m3+2 .le. d) then
c            temp1(m2+2)=x(m3+2)
c          else
c            temp1(m2+2)=0
c          end if
c3100    continue
c        call interp(ttrunc,trunc/ttrunc,temp1,dtemp1,del)
c        isign=1
c        write(*,*) 'we are'
c        call four1(trunc,isign,temp1)
c        do 3200 j=2*resterm+3,2*(trunc-resterm),1
c          temp1(j)=0
c3200    continue
c        isign=-1
c        call four1(trunc,isign,temp1)
c        do 3300 j=1,trunc,1
c          m2=2*(j-1)
c          m3=m1+(j-1)*d
c          x(m3+1)=temp1(m2+1)*invtr
c          if(m3+2.le.d) x(m3+2)=temp1(m2+2)*invtr
c3300    continue
c3000  continue
c      write(*,*) 'after',(x(i),i=1,9*d+1,d)
c**** The partial implicit method is in use. The accuracy is not ****
c**** quite essential, but the stability is very significant. Also the adaptive**
c**** scheme would be very important when the equilibrium state is close. *******
20    count1=0
      count3=0
c      nu=-0.025
      flag=0
c      do 11 while(nu.gt.-0.029.and.flag.eq.1)
      t4=5
      dh=.5E-2
      count2=0
      if(count1.eq.0) count2=-4000
c      open(3,file='torustest1.dat',status='unknown')
      do 4000 while ((t4.gt.limit.and.t4.lt.10.).and.count2.lt.2000
     \.and.dh.gt.1.e-7)
        count2=count2+1
c        if(dh.lt.0.5e-5.or.count.gt.1000) go to 10 
c        dh=min(dh,1.) 
        do 4200 i=1,tr,1
          m1=(i-1)*d
          do 4210 j=1,d,1
            a(j)=x(m1+j)
4210      continue
          call f1(d,eps1,eps2,eps3,a,dydx)  
          call jacobj(d,eps1,eps2,eps3,a,dv)        
          do 4220 j=1,d,1
            xx(m1+j)=dydx(j)
            do 4221 k=1,d,1
              jb(j,k,i)=dv(j,k)
4221        continue
4220      continue
4200    continue
c        do 4400 j=1,dtr,1
c          write(3,'(2F12.6)') x(j),xx(j)
c4400    continue  
c        stop
        t44=t4
        write(*,*) 'count2=',count2,' dh=',dh 
c        write(*,*) 'x=',(x(i),i=1,16)
        call torusfast3(omg1,omg2,x,xx,jb,y,ct2,dh,t4,
     \flag,matr,index,count3)
c        if (t44.lt.1.E-3) then
        if (t44.gt.t4*exp(0.1*dh)) then
          dh=1.01*dh
        else
          dh=2*dh/3.
c          flag=0
          t4=t44
        end if
c.and.(count3.gt.4
c     \.and.ct2.le.16)
        if(ct2.gt.16) flag=0
c        end if
        if(ct2.lt.25) then 
          do 4300 i=1,dtr,1
            x(i)=y(i)
4300      continue
        endif
c        t2=x(1)
        write(*,*) 't4=',t4,' omega=',omg1,omg2
c        write(*,*) 'count2=',count2,' dh=',dh 
4000  continue
      if(t4.le.limit) then
      count1=count1+1 
      write(*,*) 'count1=',count1
      write(4,'(3F12.6)') eps1,eps2,eps3
      do 10 i=1,dtr,1
        write(4,'(F12.6)') x(i)
10    continue
      else
        flag=0
      end if
      nu=nu-0.0001
11    continue
c      write(*,*) 'count1=',count1 
      close(4,status='keep')
      end
ccccccccccccc The following is the subroutine that will be used in CCCCCCCCCCC
ccccccccccccc the programs ccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interpfour(n,m,x01,x02,x1,x2,datk1,datk2)
      integer i,n,m,isign 
      double precision x01(n),x02(n),x1(2*n*m),x2(2*n*m)
      double precision datk1(2*n),datk2(2*n)
      double precision t1

      call dftint(n,x01,x02,datk1,datk2)
      do 1000 i=1,2*n*m,1
        x2(i)=0
1000  continue
      x2(1)=x01(1)
      x2(2)=x01(2) 
      do 2000 i=3,n/2,1
        t1=x01(i)
        x2(i)=t1
        if (mod(i,2).eq.1) then
          x2(2*n*m-i)=t1           
        else
          x2(2*n*m-i+2)=-t1
        end if
2000  continue
      isign=1
      call four1(n*m,isign,x2) 
      do 3000 i=1,n*m,1
        x1(2*i-1)=x2(2*i-1)
3000  continue
      do 4000 i=1,2*n*m,1
        x2(i)=0
4000  continue
      x2(1)=x02(1)
      x2(2)=x02(2) 
      do 5000 i=3,n/2,1
        t1=x02(i)
        x2(i)=t1
        if (mod(i,2).eq.1) then
          x2(2*n*m-i)=t1           
        else
          x2(2*n*m-i+2)=-t1
        end if
5000  continue
      isign=1
      call four1(n*m,isign,x2) 
      do 6000 i=1,n*m,1
        x1(2*i)=x2(2*i-1)
6000  continue
      end

 
      subroutine inismooth(d,trunc,x)
      integer d,trunc,i,j,m1,m2,m3,m4
      double precision t1,t2,t3,t4,h
      double precision x(d*trunc)

      t1=dble(1.)/dble(3.1415926536)
      t1=t1**5
      h=dble(2*3.1415926536)/trunc
      m1=(trunc-1)*d
      do 1000 i=1,d,1
        m2=0
        t4=h
        t2=(x(m1+i)-x(i))/2.
        t3=t2*t1
        do 1100 j=1,trunc,1
          x(m2+i)=x(m2+i)-t3*(t4-3.1416)**5
          m2=m2+d
          t4=t4+h
1100    continue
        m2=m1-d
        m3=m2-d
        m4=m3-d
        x(m3+i)=(x(m4+i)+x(m3+i)+x(m2+i))/3.
        x(d+i)=(x(i)+x(d+i)+x(2*d+i))/3.
        x(m2+i)=(x(m3+i)+x(m2+i)+x(m1+i))/3.
        x(i)=(x(m1+i)+x(i)+x(d+i))/3.
        x(m1+i)=(x(m2+i)+x(m1+i)+x(i))/3.
1000  continue
      end

      subroutine f1a(d,eps,a,dydx,t)
      integer d
      double precision eps,a(d),dydx(d)
      double precision t
      dydx(1)=a(2)
      dydx(2)=eps*(sin(a(1))+sin(a(1)-t))
      end

      subroutine jacobja(d,eps,a,btemp,dv,t)
      integer d
      double precision eps,t,t1
      double precision a(d),btemp(d,d),dv(d,d)

      t1=eps*(cos(a(1))+cos(a(1)-t))
      dv(1,1)=btemp(2,1)
      dv(1,2)=btemp(2,2)
      dv(2,1)=t1*btemp(1,1)
      dv(2,2)=t1*btemp(1,2)
      end

      subroutine f1b(d,eps,a,dydx)
      integer d
      double precision eps,a(d),dydx(d)
      dydx(2)=a(2)+eps*sin(a(1))
      dydx(1)=a(1)+dydx(2)
      end

      subroutine jacobjb(d,eps,a,dv)
      integer d
      double precision eps,a(d),dv(d,d),t
      t=eps*cos(a(1))
      dv(1,1)=1+t
      dv(1,2)=1
      dv(2,1)=t
      dv(2,2)=1
      end 

      subroutine f1(d,eps1,eps2,eps3,a,dydx)
      integer d
      double precision eps1,a(d),dydx(d)
      double precision eps2,eps3,t1
      t1=eps3*sin(a(1)+a(3))
      dydx(2)=a(2)+eps1*sin(a(1))+t1
      dydx(4)=a(4)+eps2*sin(a(3))+t1
      dydx(1)=a(1)+dydx(2)
      dydx(3)=a(3)+dydx(4)
      end

      subroutine jacobj(d,eps1,eps2,eps3,a,dv)
      integer d
      double precision eps1,eps2,eps3,t1,t2,t3
      double precision a(d),dv(d,d)

      t3=eps3*cos(a(1)+a(3))
      t1=eps1*cos(a(1))+t3
      t2=eps2*cos(a(3))+t3
      dv(1,1)=1+t1
      dv(1,2)=1
      dv(1,3)=t3
      dv(1,4)=0
      dv(2,1)=t1
      dv(2,2)=1
      dv(2,3)=t3
      dv(2,4)=0
      dv(3,1)=t3
      dv(3,2)=0
      dv(3,3)=1+t2
      dv(3,4)=1
      dv(4,1)=t3
      dv(4,2)=0
      dv(4,3)=t2
      dv(4,4)=1
      end



      subroutine interp(nn,m,dat,datemp,del)
      integer nn,m,i,j,mc
      double precision dat(2*nn),datemp(2*nn*m)
      double precision a,b,c,del,t1,t2,x2
      double precision d1,d2,d3

      do 1000 i=1,nn,1
        x2=del*(i-1)
        if (i.eq.1) then
          d1=dat(2*nn-1)
          d3=dat(3)
        elseif (i.eq.nn) then
          d1=dat(2*(nn-1)-1)
          d3=dat(1)
        else
          d1=dat(2*(i-1)-1)
          d3=dat(2*(i+1)-1)
        end if
        d2=dat(2*i-1)
        t1=(d2-d1)/del
        t2=(d3-d2)/del
        a=(t2-t1)/(2*del)
        b=t2-a*(2*x2+del)
        c=d2-a*x2**2-b*x2
        t1=del/m
        mc=(i-1)*m
        datemp(2*(mc+1)-1)=dat(2*i-1)
        do 1100 j=2,m,1
          t2=x2+(j-1)*t1
          datemp(2*(mc+j)-1)=a*t2**2+b*t2+c
1100    continue
        if (i.eq.1) then
          d1=dat(2*nn)
          d3=dat(4)
        elseif (i.eq.nn) then
          d1=dat(2*(nn-1))
          d3=dat(2)
        else
          d1=dat(2*(i-1))
          d3=dat(2*(i+1))
        end if
        d2=dat(2*i)
        t1=(d2-d1)/del
        t2=(d3-d2)/del
        a=(t2-t1)/(2*del)
        b=t2-a*(2*x2+del)
        c=d2-a*x2**2-b*x2
        t1=del/m
        datemp(2*(mc+1))=dat(2*i)
        do 1200 j=2,m,1
          t2=x2+(j-1)*t1
          datemp(2*(mc+j))=a*t2**2+b*t2+c
1200    continue
1000  continue
      do 2000 i=1,nn*m/2,1
        datemp(2*(2*i-1)-1)=datemp(2*(2*i-1)-1)*2/3
        datemp(2*(2*i-1))=datemp(2*(2*i-1))*2/3
        datemp(2*i*2-1)=datemp(2*i*2-1)*4/3
        datemp(2*i*2)=datemp(2*i*2)*4/3
2000  continue
      end






