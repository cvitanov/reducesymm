c This progam employs the first order partial implicit method to evolve
c the orbit in the phase space.  The requirement is that
c d,even, trunc,2^n  



      program  antorbfast2
      integer i,j,k,dd,trunc,tp,dtr,flag,flag1,kk,cycle,ttrunc
      parameter(dd=16,trunc=2048,dtr=dd*trunc,mm1=2*dd,mm2=2*dd)
      parameter(kk=50000,ttrunc=trunc)
      double precision tiny,hmin,eps,dh,hdid,hnext
      double precision nu,dydx(dd),yerr(dd),pts(kk*dd)
      double precision ain(dd),x(dtr),time,time1,lamold
      integer m1,m2,m3,resterm,count,ct2,ii,count1,mm3,pn
      integer index1(dtr+1),index(2*dd+1),indd(kk),indd1(kk)
      double precision yscal(dd),ytemp(dd),temp1(2*ttrunc)
      double precision xx(dtr),y(dtr),a1(dtr+1,mm1)
      double precision lam(trunc),lam1,x2(dtr+1,2*dd+1)
      double precision limit,t1,t2,t3,t4,p,invtr,L,tt(kk)
      double precision ax(dtr+1,mm1+mm2+1),aa(trunc,dd,dd)
      double precision mat(2*dd+1,2*dd+1),del,xold(dtr)
      double precision initeps,tx(ttrunc*dd),dtemp1(2*trunc)

      data tiny,hmin,eps/1.0e-30,1.0e-20,1e-4/
c  eps is used to control the tolerence error.
      pn=1
      time1=0   
      initeps=0.02
      dh=0.2e-3
c      del=2*dble(3.1415926)/ttrunc
c      nu=0.01500
      nu=1
      L=38.5
      lam1=0.134
c      lam1=4.5
      cycle=13
c      kk=40000
      resterm=2
      limit=1.e-6
      invtr=dble(1.)/trunc

     

c*************** This part is used to initialize the orbit***************
c***Fourier transform the evolution orbit then erase the high-f part **** 
500   open(2,file='antorb4b.dat',status='old')
      open(4,file='antcyc4e3.dat',status='unknown')      
c      read(2,'(F12.8,2F13.8)') ain(1),ain(2),ain(3)
c      write(*,'(F12.8,2F13.8)') ain(1),ain(2),ain(3)
c      stop
      i=0
      do 1000 while(i.lt.50000*dd)
        i=i+1
        read(2,'(F12.6)') pts(i)
1000  continue
c      do 1500 i=1,100,1
c        write(*,*) pts(1+(100*i+25000)*dd)
c1500  continue
      close(2,status='keep')
c      go to 20
c      stop


      count=0
      i=19000
      m3=24000
      do 1600 while(i.lt.21000)
        i=i+2
        cycle=21500-i
        do 1610 while (cycle .lt. min(2000,m3-i))
          cycle=cycle+1
          t2=0
          t3=0
          m1=(i-1)*dd
          m2=(i+cycle-1)*dd
          do 1611 j=1,dd,1
            t2=t2+(pts(m2+j)+pts(m2+j+dd)-pts(m1+j)
     \-pts(m1+j+dd))**2
            t3=t3+pts(m2+j)**2+pts(m1+j)**2
c            t2=t2+((pts(m2+dd+j)-pts(m2-dd+j)-pts(m1+dd+j)
c     \+pts(m1-dd+j))/0.2)**2
1611      continue
          t2=dsqrt(t2/(2*t3))
          if (t2.lt.initeps .and. count.lt.kk) then
            count=count+1
            indd(count)=i
            indd1(count)=cycle
            cycle=cycle+15
            i=i+1
c          write(*,*) 'The ',count,'th one:index=',i
c          write(4,'(3F12.6)') pts(count,1),pts(count,2),pts(count,3)
c          write(4,'(2F12.6)') pts(count,4),pts(count,5)         
          endif
1610    continue
1600  continue
c      close(4,status='keep')
c      write(*,*) 'count=',count
      if (count.eq.0) then
        write(*,*) 'There is no count.'
        stop
      endif
c      stop
CCCCCFinish matching time series and start to initialize the dataCCCCC

20    count1=0
      mm3=count 
      write(*,*) 'count=',count 
      do 10,ii=1,mm3,max(1,mm3/200)
      m1=(indd(ii)-1)*dd
c       m1=1580*dd 
      time=indd1(ii)*0.2
c      time=12.
      write(*,*) 'ii=',ii,' m3=',mm3
      write(*,*) 'time=',time,' count1=',count1 

      time1=time/trunc
      lam1=time/(2*3.1415927)
      resterm=time/4
c      resterm=16
c      m1=(13592-1)*dd
c      do 11 i=1,indd1(ii)+1,1
c        do 12 j=1,dd,1
c          write(4,'(F12.6)') pts(m1+j)
c12      continue
c        m1=m1+dd
c11    continue
c      close(4,status='keep')
c      stop
      do 1800 i=1,dd,1
        ain(i)=pts(m1+i)
1800  continue
c      tp=1
c      p=0.22
c      call antprop(L,dd,nu,ain,dydx,yscal,yerr,ytemp,18*time,tp,p,pn)
      do 2000 i=1,trunc,1
c        write(*,*) 'i=',ain
        m1=(i-1)*dd
        tp=1
        p=0.22
c         call propagat(dd,nu,ain,dydx,yscal,yerr,ytemp,time1,tp,p)
        call antprop(L,dd,nu,ain,dydx,yscal,yerr,ytemp,time1,tp,p,pn)
        do 2100 j=1,dd,1
          x(m1+j)=ain(j)
2100    continue
2000  continue
      count2=2
c      go to 30

40    write(*,*) 'before',(x(i),i=1,9*dd+1,dd)
c      call inismooth(dd,trunc,x)
c      if (count2.eq.5) then
c        write(*,*) ' done without findings'
c        go to 10
c      end if
c      count2=count2+1
      do 3000 i=1,dd/2,1
          m1=2*(i-1)
        do 3100 j=1,trunc,1
          m2=2*(j-1)
          m3=m1+(j-1)*dd
          temp1(m2+1)=x(m3+1)
          temp1(m2+2)=x(m3+2)
3100    continue
c        call interp(ttrunc,trunc/ttrunc,temp1,dtemp1,del)
        isign=1
c        write(*,*) 'we are'
        call four1(trunc,isign,temp1)
c        write(*,*) 'we are 2',dtemp1(1),dtemp1(2),dtemp1(3)
        do 3200 j=2*resterm+3,2*(trunc-resterm),1
          temp1(j)=0
3200    continue
c        if (i.eq.mod(count2,9)) then
c          temp1(trunc+1)=(0.0015*i*(6-i)+0.03)/(count2-2)**3
c          temp1(trunc+2)=(0.002*i*(6-i)+0.02)/(count2-2)**3
c        end if
        isign=-1
        call four1(trunc,isign,temp1)
        do 3300 j=1,trunc,1
          m2=2*(j-1)
          m3=m1+(j-1)*dd
          x(m3+1)=temp1(m2+1)*invtr
          x(m3+2)=temp1(m2+2)*invtr
3300    continue
3000  continue
      write(*,*) 'after',(x(i),i=1,9*dd*trunc/ttrunc+1,dd*trunc/ttrunc)
c      resterm=30
c      stop
c**** The partial implicit method is in use. The accuracy is not ****
c**** quite essential, but the stability is very significant. Also the adaptive**
c**** scheme would be very important when the equilibrium state is close. *******
30    t1=50
      flag=1
      flag1=0 
      count=0
      count2=0
      ct2=100 
      dh=1.0e-3
      t4=1.05*t1
      lamold=lam1
      do 3900 i=1,dtr,1
        xold(i)=x(i)
3900  continue
      do 4000 while (t4.gt.limit.and.t1.lt.100)
c        count2=count2+1
c        t4=t1
        if(dh.lt.0.5e-5.or.count.gt.1000) go to 10 
        dh=min(dh,1.) 
        call antfast(L,nu,dh,x,xx,y,lam,flag,aa,ax,a1,
     \count,lam1,t1,ct2,index1,index,x2,mat,flag1)
c        ct2=20
        if (ct2.gt.16) then
          t2=flag1-6.5+(16-ct2)/8.0
c          if (t2.ge.0) dh=dh*(1+t2*0.08)
          if (t2.lt.0) dh=dh*(1+t2**3*0.0033)
        endif
        if (t4.gt.t1) then
          do 4100 i=1,dtr,1
            xold(i)=x(i)
            x(i)=y(i)
c            write(3,'(F12.6)') y(i)
4100      continue
          t4=1.05*t1
          lamold=lam1
          count=count+1
          if(ct2.gt.3) dh=max(2.e-4,dh)
        else
          do 4200 i=1,dtr,1
            x(i)=xold(i)
4200      continue
          lam1=lamold
          dh=0.3*dh
          flag1=5
          ct2=18
c          go to 40
        end if
        if (flag1.gt.8) then
          t2=flag1-8
          dh=dh*(1+0.2)
        end if


c        t1=0
c        do 4200 i=1,trunc,1
c          t1=dmax1(t1,lam(i))
c4200    continue
        write(*,*) 'count=' ,count,' iter=',ct2,'dh=',dh
        write(*,*) 'lambda, i=1,10',(lam(i),i=1,10)
        write(*,*) 'x value',(x(i),i=1,9*dd+1,dd)
4000  continue
      if (t1.le.limit) then
        i=0
        flag=0
        do 4550 while(i.lt.count1.and.flag.eq.0)
          i=i+1
          if(abs(lam1-tt(i))/lam1.lt.1e-4) flag=1
4550    continue
        if (flag.eq.0) then 
          write(4,'(F18.12)') lam1
          do 4500 i=1,dtr,1
            write(4,'(F18.12)') x(i)
4500      continue
          count1=count1+1
          tt(count1)=lam1
        end if
      end if
10    continue
      write(*,*) 'count1=',count1 
      close(4,status='keep')
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









