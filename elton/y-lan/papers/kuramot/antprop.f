c This program is used to propagate the point a(trunc) for the map on a1=0
c Poincare section or the corresponding point for the Jacobian according to 
c the value of the argument tp{0,1}. The modified adaptive fifth-order  
c Runge-Kutta method (runge)is used  Please pay attention to the arbitrarily 
c chosen dimensions for the transient array aks below.  



      subroutine  antprop(L,trunc,nu,a,dydx,yscal,yerr,ytemp,time,
     \tp,p,pn)
      integer trunc,tp,pn
      double precision tiny,hmin,eps,htry,hdid,hnext,p
      double precision nu,dydx(trunc),yerr(trunc),L
      double precision a(trunc),b(trunc),time,time1
      integer i,count
      double precision yscal(trunc),ytemp(trunc)
      data tiny,hmin,eps/1.0e-30,1.0e-29,1.0e-6/
c  eps is used to control the tolerence error.
      time1=0
      htry=1.0e-4
      count=0

500   call antks(L,trunc,nu,a,dydx)
c      count=count+1
c      if(count.eq.12) stop
c      write(*,*) 'dydx',dydx(1),dydx(6),dydx(20)
c      write(*,*) 'v of a',a(1),a(6),a(20)
c      write(*,*) 'htry',htry
      do 1000 i=1,trunc
        yscal(i)=dabs(a(i))+dabs(dydx(i)*htry)+tiny
1000  continue
      call rkqs(L,trunc,nu,a,time1,dydx,yscal,yerr,
     \ytemp,eps,htry,hdid,hnext,tp,p,pn)
      if(time1.gt.350) return
      if (tp.eq.3) then
        time=time1
c        write(*,*) 'hdid',hdid
        tp=0
        return
      end if
      if(tp.eq.1.and.time1+hnext.gt.time) hnext=time-time1
      if(tp.eq.1.and.time1.ge.time) return
      if (hnext.le.hmin) then
        write(*,*) 'step is too small in propagat'
        stop
      end if
      htry=hnext
      go to 500
      
      end




c This routine needs the initial value of x,y. n is the dimension of out
c equation set.yscal,yerr and eps are error control variables.
      subroutine rkqs(L,n,nu,y,x,dydx,yscal,yerr,ytemp,eps
     \,htry,hdid,hnext,tp,p,pn)
      integer n,tp,pn
      double precision y(n),dydx(n),nu,L
      double precision yscal(n),yerr(n),ytemp(n)
      double precision eps,htry,hdid,hnext
      double precision h,htemp,p,temp,delt
      real safty,pgrow,pshk,errcon
      double precision x,xnew,errmax
      integer i,flag,count
      safty=0.9
      pgrow=-0.2
      pshk=-0.25
      errcon=1.89e-4
      count=0
c here is the relation errcon=(5/safty)**(1/pgrow)

      flag=0
      h=htry
      do 1000 while(flag.eq.0)
        call rkck(L,n,nu,y,x,dydx,h,yerr,ytemp)
        count=count+1
        errmax=0
        do 1100 i=1,n
          errmax=dmax1(errmax,dabs(yerr(i)/yscal(i)))
1100    continue
c        write(*,*) 'errmax',errmax,yscal(1),yscal(3)
        errmax=errmax/eps
        if (errmax.le.1) then 
          if (y(pn).lt.p.and.ytemp(pn).gt.p.and.tp.eq.0) then
            temp=p-y(pn)
            call rkck2(L,n,nu,y,x,dydx,temp,yerr,ytemp,delt,pn) 
            hdid=delt
            x=x+hdid
            hnext=h
            do 1200 i=1,n
              y(i)=ytemp(i)
1200        continue
            y(pn)=p 
            tp=3
            go to 3000 
          else  
            flag=1
            go to 1000
          end if  
        end  if 
        htemp=safty*h*(errmax**pshk)
        if(htemp.lt.0.1*h) htemp=0.1*h
        h=htemp
        xnew=x+h
        if (xnew.eq.x) then
          write(*,*) 'stepsize underflow in rkqs'
          stop
        endif
1000  continue
      if (errmax.gt.errcon) then 
        hnext=safty*h*(errmax**pgrow)     
      else
        hnext=5.0*h
      end if
c      write(*,*) 'h,hnext,errmax,hnext',h,hnext,errmax,hnext
c      write(*,*) 'count',count
c      stop
      hnext=min(2e-2,hnext)
      hdid=h
      x=x+h
      do 2000 i=1,n
        y(i)=ytemp(i)
2000  continue
3000  continue 
      end   

      subroutine rkck(L,n,nu,y,x,dydx,h,yerr,yout) 
      integer n,trunc
c Here comes the awkward part of the subroutine since you have to give a large 
c dimension to the transient array. It is 50 here.
      parameter (trunc=16)
      double precision y(n),yout(n),L
      double precision x,h,dydx(n),nu
      double precision ytemp(trunc),yerr(n)
      double precision a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43
      double precision b51,b52,b53,b54,b61,b62,b63,b64,b65
      double precision c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6
      double precision ak2(trunc),ak4(trunc),ak3(trunc),ak5(trunc)
      double precision ak6(trunc),ak1(trunc)
      integer i
      data a2,a3,a4,a5,a6,b21/0.2,0.3,0.6,1.0,0.875,0.2/
      data b31,b32,b41,b42,b43/0.075,0.225,0.3,-0.9,1.2/
      data b51,b52,b53,b54,b61,b62,b63,b64,b65/-0.2037037037037,2.5,
     \-2.5925925925926,1.296296296296,0.0294958043981,0.341796875,
     \0.0415943287037,0.4003454137731,0.061767578125/
      data c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6/0.0978835978836,
     \0.402576489533,0.2104377104377,0.2891022021457,-0.0042937748016,
     \.0186685860938,-0.0341550268308,-0.0193219866071,.0391022021457/
      
      
      
      do 900 i=1,n
        ak1(i)=dydx(i)
900   continue
      do 1000 i=1,n
        ytemp(i)=y(i)+b21*h*ak1(i)
1000  continue
      call antks(L,n,nu,ytemp,ak2)
      do 2000 i=1,n
        ytemp(i)=y(i)+h*(b31*ak1(i)+b32*ak2(i))
2000  continue
      call antks(L,n,nu,ytemp,ak3)
      do 3000 i=1,n
        ytemp(i)=y(i)+h*(b41*ak1(i)+b42*ak2(i)+b43*ak3(i))
3000  continue
      call antks(L,n,nu,ytemp,ak4)
      do 4000 i=1,n
        ytemp(i)=y(i)+h*(b51*ak1(i)+b52*ak2(i)+b53*ak3(i)+b54*ak4(i))
4000  continue
      call antks(L,n,nu,ytemp,ak5)
      do 5000 i=1,n
        ytemp(i)=y(i)+h*(b61*ak1(i)+b62*ak2(i)+b63*ak3(i)+
     \b64*ak4(i)+b65*ak5(i))
5000  continue
      call antks(L,n,nu,ytemp,ak6)
      do 6000 i=1,n
        yout(i)=y(i)+h*(c1*ak1(i)+c3*ak3(i)+c4*ak4(i)+c6*ak6(i))
        yerr(i)=h*(dc1*ak1(i)+dc3*ak3(i)+dc4*ak4(i)
     \+dc5*ak5(i)+dc6*ak6(i))
6000  continue
      end
      
      subroutine rkck2(L,n,nu,y,x,dydx,h,yerr,yout,time,pn) 
      integer n,trunc,pn
c Here comes the awkward part of the subroutine since you have to give a large 
c dimension to the transient array. It is 50 here.
      parameter (trunc=16)
      double precision y(n),yout(n),L
      double precision x,h,dydx(n),nu,time
      double precision ytemp(trunc),yerr(n)
      double precision a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43
      double precision b51,b52,b53,b54,b61,b62,b63,b64,b65
      double precision c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6
      double precision ak2(trunc),ak4(trunc),ak3(trunc),ak5(trunc)
      double precision ak6(trunc),ak1(trunc)
      double precision t1,t2,t3,t4,t5,t6
      integer i
      data a2,a3,a4,a5,a6,b21/0.2,0.3,0.6,1.0,0.875,0.2/
      data b31,b32,b41,b42,b43/0.075,0.225,0.3,-0.9,1.2/
      data b51,b52,b53,b54,b61,b62,b63,b64,b65/-0.2037037037037,2.5,
     \-2.5925925925926,1.296296296296,0.0294958043981,0.341796875,
     \0.0415943287037,0.4003454137731,0.061767578125/
      data c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6/0.0978835978836,
     \0.402576489533,0.2104377104377,0.2891022021457,-0.0042937748016,
     \.0186685860938,-0.0341550268308,-0.0193219866071,.0391022021457/
      
      
      t1=1/dydx(pn)
c      ak1(1)=1
      do 900 i=1,n
        ak1(i)=dydx(i)*t1
900   continue
      ak1(pn)=1
      do 1000 i=1,n
        ytemp(i)=y(i)+b21*h*ak1(i)
1000  continue
      call antks(L,n,nu,ytemp,ak2)
      t2=1/ak2(pn)
c      ak2(1)=1
      do 1500 i=1,n
        ak2(i)=ak2(i)*t2
1500  continue
      ak2(pn)=1
      do 2000 i=1,n
        ytemp(i)=y(i)+h*(b31*ak1(i)+b32*ak2(i))
2000  continue
      call antks(L,n,nu,ytemp,ak3)
      t3=1/ak3(pn)
c      ak3(1)=1
      do 2500 i=1,n
        ak3(i)=ak3(i)*t3
2500  continue
      ak3(pn)=1
      do 3000 i=1,n
        ytemp(i)=y(i)+h*(b41*ak1(i)+b42*ak2(i)+b43*ak3(i))
3000  continue
      call antks(L,n,nu,ytemp,ak4)
      t4=1/ak4(pn)
c      ak4(1)=1
      do 3500 i=1,n
        ak4(i)=ak4(i)*t4
3500  continue
      ak4(pn)=1
      do 4000 i=1,n
        ytemp(i)=y(i)+h*(b51*ak1(i)+b52*ak2(i)+b53*ak3(i)+b54*ak4(i))
4000  continue
      call antks(L,n,nu,ytemp,ak5)
      t5=1/ak5(pn)
c      ak5(1)=1
      do 4500 i=1,n
        ak5(i)=ak5(i)*t5
4500  continue
      ak5(pn)=1
      do 5000 i=1,n
        ytemp(i)=y(i)+h*(b61*ak1(i)+b62*ak2(i)+b63*ak3(i)+
     \b64*ak4(i)+b65*ak5(i))
5000  continue
      call antks(L,n,nu,ytemp,ak6)
      t6=1/ak6(pn)
c      ak6(1)=1
      do 5500 i=1,n
        ak6(i)=ak6(i)*t6
5500  continue
      ak6(pn)=1
      do 6000 i=1,n
        yout(i)=y(i)+h*(c1*ak1(i)+c3*ak3(i)+c4*ak4(i)+c6*ak6(i))
        yerr(i)=h*(dc1*ak1(i)+dc3*ak3(i)+dc4*ak4(i)
     \+dc5*ak5(i)+dc6*ak6(i))
6000  continue
      time=h*(c1*t1+c3*t3+c4*t4+c6*t6)
      end   
     









