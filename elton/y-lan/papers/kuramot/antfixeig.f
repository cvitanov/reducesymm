c This program is used to give the Poincare sections which pass the
c fixed points and cut the plane of unstable directions.
 
      program antfixeig
      integer d,trunc,i,j,k,n,m1,m
      parameter (n=4,d=16,trunc=512)
      integer isign,index(d),check
      double complex mat(d,d),av(d),vv(d),af(d,n)
      double precision orb(trunc),dv(d,d),temp1(2*trunc),a(d)
      double precision orbdt(n*trunc,4),eir(d),eii(d),vv1(d)
      double precision L,nu,pi,t1,t2,t3,b(d)
      double precision ct,st,h,ap(d,n)
      double precision xold(d),g(d),p(d),fjac(d,d),fvec(d)

      L=38.5
      nu=1
      h=L/trunc
      pi=dble(3.1415926536)
      open(1,file='antcycl7.dat',status='old')
c      open(2,file='antfix4b1.dat',status='old')
c      open(1,file='antorb4b.dat',status='old')
      open(2,file='antfix4b12.dat',status='old')
      do 500 i=1,d*2500,1
        read(1,'(F12.6)') t1
500   continue
      do 600 i=1,d,1
        read(1,'(F12.6)') b(i)
600   continue
      do 1000 i=2,2,1
        m1=0
        do 1100 j=1,trunc*n,1
          m1=m1+1
          read(i,'(F12.6)') orbdt(m1,i)
1100    continue
1000  continue
      do 1500 i=1,trunc*n*6,1
        read(2,'(F12.6)') t1
1500  continue
      do 1600 i=1,2,1
        do 1610 j=1,trunc*n,1
          read(2,'(F12.6)') orbdt(j,i+2)
1610    continue
1600  continue
      close(1,status='keep')
      close(2,status='keep')
c      close(3,status='keep')
c      close(4,status='keep')
c      write(*,*) '1st place'  
      open(1,file='antfix4p1.dat',status='unknown')
      open(2,file='antfix4p11.dat',status='unknown')
c      open(3,file='antfix4p12.dat',status='unknown')
c      open(4,file='antfix4p13.dat',status='unknown')
c      write(*,*) '2nd place'      
      do 2000 i=1,4,1
        write(*,*) 'i=',i
        if (i.ne.1) then
        m1=0
        do 2100 j=1,trunc*n,n
          m1=m1+1
          orb(m1)=orbdt(j,i)
          temp1(2*m1-1)=orb(m1)
          temp1(2*m1)=0
2100    continue
        isign=-1
        call four1(trunc,isign,temp1)
        t1=atan2(temp1(4),temp1(3))-pi/2
c        write(*,*) '2nd place'
c         t1=
        do 2200 j=1,d,1
          m1=2*j+1       
          ct=cos(j*t1)
          st=sin(j*t1)
          t2=temp1(m1+1)*ct-temp1(m1)*st
          a(j)=t2/trunc 
c          b(j)=(temp1(m1+1)*st+temp1(m1)*ct)/trunc
2200    continue
        else
        do 2250 j=1,d,1
          a(j)=b(j)
2250    continue
        end if
        call newt(a,xold,g,p,fjac,fvec,d,nu,L,check,index,vv1,t1)
c        write(*,*) 'b=',b
        call antks(L,d,nu,a,vv1)
        call ksjcb(L,d,nu,dv,a)
        write(*,*) '2nd place=',a
        do 2300 j=1,d,1
          ap(j,i)=a(j) 
          do 2310 k=1,d,1
            mat(j,k)=cmplx(dv(j,k),0) 
2310      continue
2300    continue       
        call balanc(dv,d)
        call hessen(dv,d)
        call hqr(dv,d,eir,eii)
        write(*,*)'The result',eir
        t1=-1000.
        do 2400 j=1,d,1
          if (t1.lt.eir(j)) then
            t1=eir(j)
            m=j
          end if
2400    continue
        do 2500 j=1,d,1
          mat(j,j)=mat(j,j)-cmplx(eir(m),eii(m))
2500    continue
        call heigc(mat,d,av,index,vv,vv1)
        do 2600 j=1,d,1
          af(j,i)=av(j)
2600    continue
2000  continue
      t2=0
      t3=0
      do 3000 i=1,d,1
        t1=ap(i,1)-ap(i,2)
        t2=t2+t1*dble(af(i,2))
        t3=t3+t1*imag(af(i,2))
3000  continue
      t1=(0-t2)/t3
      t3=0
      do 4000 i=1,d,1
        t2=dble(af(i,2))+t1*imag(af(i,2))
        t3=t3+ap(i,2)*t2
        write(1,'(F12.6)') t2 
4000  continue
      write(1,'(F12.6)') t3
      t2=0
      t3=0
      do 5000 i=1,d,1
        t1=ap(i,4)-ap(i,3)
        t2=t2+t1*dble(af(i,3))
        t3=t3+t1*imag(af(i,3))
5000  continue
      t1=(0-t2)/t3
      t3=0
      do 6000 i=1,d,1
        t2=dble(af(i,3))+t1*imag(af(i,3))
        t3=t3+ap(i,3)*t2
        write(2,'(F12.6)') t2 
6000  continue
      write(2,'(F12.6)') t3
      close(1,status='keep')
      close(2,status='keep')
c      close(3,status='keep')
c      close(4,status='keep')
      write(*,*) 'done'
      end

CCCCCCCCCCCCCCCCCCCCCCCC Root-search Routine belowCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine newt(x,xold,g,p,fjac,fvec,n,nu,L,check,index,vv,f)
      integer n,check
      double precision x(n),g(n),p(n),xold(n)
      double precision fjac(n,n),vv(n),temp1,L
      double precision fvec(n),tolf,tolmin,tolx,stpmx,nu
      double precision d,den,f,fold,stpmax,sum,temp,test,func
      integer i,j,its,nn,index(n),maxits
      data maxits,tolf,tolmin,tolx,stpmx/150,1.e-5,1.e-6,1.e-7,100./
      nn=n
      f=func(nn,nu,L,x,fvec)
      test=0
      do 1000 i=1,n
        temp=dabs(fvec(i))
        if(temp.gt.test) test=temp
1000  continue
      if (test.lt.0.01*tolf) then
        check=0
        RETURN
      end if
      sum=0
      do 2000 i=1,n
        sum=sum+x(i)*x(i)
2000  continue
      stpmax=stpmx*dmax1(dsqrt(sum),dble(n))
      do 3000 its=1,maxits
        call antks(L,n,nu,x,fvec)
        call ksjcb(L,n,nu,fjac,x)
c        call f1(n,nu,x,fvec)
c        call jacob(n,nu,x,fjac)
c        write(*,*) 'fjac and fvec',fjac,fvec 
        do 3100 i=1,n
          sum=0
          do 3110 j=1,n
            sum=sum+fjac(j,i)*fvec(j)
3110      continue
          g(i)=sum
          p(i)=-fvec(i)
          xold(i)=x(i)
3100    continue
        fold=f   
        call lucmp(n,fjac,index,vv)
        call lusub(n,fjac,p,index)
c        write(*,*) 'the previous',p 
        call lintrk(n,xold,fold,g,p,x,f,fvec,stpmax,check,nu,L)
c        write(*,*) 'the value of f',f
        test=0
        do 3200 i=1,n 
          temp=dabs(fvec(i))
          if(test.lt.temp) test=temp
3200    continue
        if (test.lt.tolf) then
          check=20
          RETURN
        end if
        if (check.eq.1) then
          test=0
          den=dmax1(f,dble(0.5*n))
          temp1=0
          do 3300 i=1,n
            temp=dabs(g(i))*dmax1(dabs(x(i)),dble(1.0))/den
            if(temp.gt.test) test=temp
            temp=dabs(fvec(i))
            if(temp.gt.temp1) temp1=temp
3300      continue
          if (test.lt.tolmin.or.temp1.gt.tolf) then
            check=1
          else
c            write(*,*)'g(i)',g,x,den,test,tolmin
            check=30
          end if   
          RETURN
        end if
        test=0
        do 3400 i=1,n
          temp=dabs(x(i)-xold(i))/dmax1(dabs(x(i)),dble(1.0))
          if(temp.gt.test) test=temp
3400    continue  
        if (test.lt.tolx) then 
           check=40        
c          write(*,*)'x value',x
c          write(*,*)'xold value',xold
          RETURN
        end if
3000  continue
c      write(*,*)'maxits exceeded in newt'
      end  


      subroutine lintrk(n,xold,fold,g,p,x,f,fvec,stpmax,check,nu,L)
      integer n,check
      double precision xold(n),fold,g(n),p(n),x(n),L
      double precision stpmax,alf,tolx,func,f,nu,fvec(n) 
      double precision a,alam,alam2,alamin,b,disc,f2,fold2,alf1
      double precision rhs1,rhs2,slope,sum,temp,test,tmplam,fbig
      integer i,j
      data alf,alf1,tolx,fbig /1.0e-4,1e-2,1.0e-7,.5/
     
      check=0
      sum=0
      do 1000 i=1,n
        sum=sum+p(i)*p(i)
1000  continue 
      sum=dsqrt(sum)
      if (sum.gt.stpmax) then
        do 2000 i=1,n
          p(i)=p(i)*stpmax/sum
2000    continue
      end if
      slope=0 
      do 3000 i=1,n
        slope=slope+g(i)*p(i)
3000  continue
      test=0
      do 4000 i=1,n
        temp=dabs(p(i))/dmax1(dabs(xold(i)),dble(1.0))
        if(temp.gt.test) test=temp
4000  continue 
      alamin=tolx/test
      alam=1.0
      do 5000 while(check.eq.0)
        if (alam.lt.alf1.and.f.gt.fbig) then
          do 5500 i=1,n
            x(i)=xold(i)-alam*g(i)
5500      continue
        else 
          do 5100 i=1,n
            x(i)=xold(i)+alam*p(i)
5100      continue
        end if
        f=func(n,nu,L,x,fvec)
c        write(*,*)'f and fold',f,fold
        if (alam.lt.alamin) then
          do 5200 i=1,n
            x(i)=xold(i)
5200      continue
          check=1          
          RETURN          
        else if (f.le.fold+alf*alam*slope) then          
          RETURN
        else
          if(alam.eq.1.0) then
            tmplam=-slope/(2.0*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope 
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2)
            b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/
     \(alam-alam2)
            if (a.eq.0) then 
              tmplam=-slope/(2.0*b)
            else  
              disc=b*b-3.0*a*slope
              if (disc.lt.0) then
                write(*,*) 'roundoff error in lintrk'
              else
                tmplam=(-b+dsqrt(disc))/(3.0*a)
              end if
            end if
            if(tmplam.gt.0.5*alam) tmplam=0.5*alam 
          end if
        end if
        alam2=alam
c        write(*,*) 'alam',alam,x(1),slope,p(1),p(20),p(10)
c        write(*,*) 'p and g',p,g
        f2=f
        fold2=fold
        alam=dmax1(tmplam,0.1*alam)
5000  continue
      end      

      function func(nn,nu,L,x,fvec)
      integer i,nn
      double precision nu,L,x(nn),fvec(nn)
      double precision func

      call antks(L,nn,nu,x,fvec) 
      func=0
      do 1000 i=1,nn,1
        func=func+fvec(i)**2
1000  continue
      func=func*0.5
      end
