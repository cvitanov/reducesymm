c This routine is using variational method to find 2-d tori in a general flow 
c provided that the initial guess is quite close to a true torus.

      subroutine torusfast(omg,x,xx,jb,y,ct2,dh,t8,t88)
      integer d,trunc,dtr,m1,m2,m3,m4,m5,m22,m33,ct2,ktrunc,kdtr
      parameter(d=16,trunc=128,dtr=d*trunc,ktrunc=trunc/2)
      parameter(kdtr=dtr/2)
      integer i,j,k,m,n,isign,index(kdtr+1),indt(2)
      double precision omg,x(dtr),xx(dtr),jb(d,d,trunc),y(dtr)
      double precision jbfr(ktrunc,d,d),matr(kdtr+1,kdtr+1)
      double precision tempk1(trunc),tempk2(trunc)
      double precision t1,t2,t3,t4,t5,t8,t88,dh,pi,t0
      double precision a(kdtr),b(kdtr),dx(kdtr+1)
      double precision datk1(2*trunc),datk2(2*trunc)
      double precision dx0(kdtr+1),vt1(kdtr+1),matr0(kdtr+1,kdtr+1)
     
      pi=3.1415926535897932D0
      t0=2*pi/trunc
      n=ktrunc/2-1
      do 500 i=1,kdtr+1,1
        dx(i)=0
        dx0(i)=0
        do 510 j=1,kdtr+1,1
          matr(i,j)=0
510     continue
c        if(mod(i,2).eq.0.and.i.le.10) matr(kdtr+1,i)=1
500   continue

cccccc From x,xx,jb, generate the Fourier coefficients a,b,jbfr ccccccccccccc
      t88=0
      do 1000 i=1,d,1
        do 1100 j=1,trunc,1
          m1=(j-1)*d+i
          t1=(j-1)*t0
          if (mod(i,2).eq.-10) then
            tempk1(j)=x(m1)-t1
            tempk2(j)=xx(m1)-t1
          else
            tempk1(j)=x(m1)
            tempk2(j)=xx(m1)
          end if
1100    continue
        call dftint(trunc,tempk1,tempk2,datk1,datk2)
        do 1200 j=1,ktrunc,1
          m1=(j-1)*d+i
          a(m1)=tempk1(j)
          b(m1)=tempk2(j)
c          t88=t88+(a(m1)-b(m1))**2
1200    continue
1000  continue
c      write(*,*)'a=',(a(j),j=1,10),' b=',(b(j),j=1,10)
c      write(*,*)'x=',(x(j),j=1,10),' xx=',(xx(j),j=1,10)
      if (mod(d,2).eq.1) then
        m1=d*d-1
        do 2000 i=1,trunc,1  
          tempk1(i)=jb(d,d,i)
          tempk2(i)=0
2000    continue
        call dftint(trunc,tempk1,tempk2,datk1,datk2)
        do 2500 i=1,ktrunc,1
          jbfr(i,d,d)=tempk1(i)
2500    continue
      else 
         m1=d*d
      end if
      do 3000 i=1,m1,2
        m2=1+(i-1)/d
        m3=mod(i,d)
        if(m3.eq.0) m3=d
        m4=1+i/d
        m5=mod(i+1,d)
        if(m5.eq.0) m5=d   
        do 3100 j=1,trunc,1
          tempk1(j)=jb(m2,m3,j)
          tempk2(j)=jb(m4,m5,j) 
3100    continue
        call dftint(trunc,tempk1,tempk2,datk1,datk2)
        do 3200 j=1,ktrunc,1
          jbfr(j,m2,m3)=tempk1(j)
          jbfr(j,m4,m5)=tempk2(j)
3200    continue
3000  continue

cccccc Construct the big matrix matr in the variational equation ccccccccccc
       t8=0
       do 4000 k=0,n,1
         m1=2*k*d
         t1=cos(k*omg)
         t2=sin(k*omg)
         do 4100 j=1,d,1
           m2=m1+j
           m3=m1+d+j
           matr(m2,m2)=t1
           matr(m2,m3)=-t2
           matr(m2,kdtr+1)=-k*(a(m2)*t2+a(m3)*t1)
c           if(k.eq.0.and.(j.eq.1.or.j.eq.3)) matr(m2,kdtr+1)
c     \= matr(m2,kdtr+1)+1 
           matr(m3,m2)=t2
           matr(m3,m3)=t1
           matr(m3,kdtr+1)=k*(a(m2)*t1-a(m3)*t2)
c           if(j.eq.4) then
c           if(k.eq.0) then
c               matr(kdtr+1,m2)=1
c           else
c               matr(kdtr+1,m2)=2*cos(k*6*t0)
c               matr(kdtr+1,m3)=-2*sin(k*6*t0)
c           end if
c           end if
           matr(kdtr+1,m2)=-k*a(m3)
           matr(kdtr+1,m3)=k*a(m2)
c           matr(kdtr+1,m2)=a(m2)
c           matr(kdtr+1,m3)=a(m3)
c           if (k.eq.0) then
c             matr(kdtr+1,m2)=0.5*a(m2)
c             matr(kdtr+1,m3)=0.5*a(m3)
c           end if
4100     continue         
         do 4200 m=1,n-k,1
           do 4210 j=1,d,1
             m4=m1+j
             m5=m4+d 
             m2=2*(m+k)*d
             m22=2*m*d
             m3=m2+d
             m33=m22+d
             do 4211 i=1,d,1
               matr(m4,m22+i)=matr(m4,m22+i)-jbfr(2*(m+k)+1,j,i)
               matr(m4,m33+i)=matr(m4,m33+i)-jbfr(2*(m+k)+2,j,i)
               matr(m4,m2+i)=matr(m4,m2+i)-jbfr(2*m+1,j,i)
               matr(m4,m3+i)=matr(m4,m3+i)-jbfr(2*m+2,j,i)
               matr(m5,m22+i)=matr(m5,m22+i)-jbfr(2*(m+k)+2,j,i)
               matr(m5,m33+i)=matr(m5,m33+i)+jbfr(2*(m+k)+1,j,i)
               matr(m5,m2+i)=matr(m5,m2+i)+jbfr(2*m+2,j,i)
               matr(m5,m3+i)=matr(m5,m3+i)-jbfr(2*m+1,j,i)
4211         continue             
4210       continue
4200     continue
         do 4300 m=0,k,1
           do 4310 j=1,d,1
             m4=m1+j
             m5=m4+d
             m22=2*m*d
             m33=m22+d
             do 4311 i=1,d,1
               matr(m4,m22+i)=matr(m4,m22+i)-jbfr(2*(k-m)+1,j,i)
               matr(m4,m33+i)=matr(m4,m33+i)+jbfr(2*(k-m)+2,j,i)
               matr(m5,m22+i)=matr(m5,m22+i)-jbfr(2*(k-m)+2,j,i)
               matr(m5,m33+i)=matr(m5,m33+i)-jbfr(2*(k-m)+1,j,i)
4311         continue
4310       continue
4300     continue
         do 4400 j=1,d,1
           m2=m1+j
           m3=m2+d
           t3=b(m2)-(a(m2)*t1-a(m3)*t2)
c           if(k.eq.0.and.(j.eq.1.or.j.eq.3)) t3=t3-omg 
           t4=b(m3)-(a(m2)*t2+a(m3)*t1)
           dx(m2)=t3
           dx(m3)=t4
           dx0(m2)=t3
           dx0(m3)=t4
           t88=t88+t3**2+t4**2
           t5=max(abs(t3),abs(t4))
           if (t5.gt.t8) then         
             t8=t5  
c             indt(1)=k
c             indt(2)=j
           end if
c           write(*,*) 't3,t4=',t3,t4                    
4400     continue
4000   continue
c       write(*,*) 't8=',t8, ' ind=',indt
c       stop 
c       matr(kdtr+1,kdtr+1)=1
c       matr(kdtr+1,2*d+1)=1

cccccc Solve the linear equation for the change of the Fourier coeffient cccc
       do 5000 i=1,kdtr+1,1
         do 5100 j=1,kdtr+1,1
           matr0(i,j)=matr(i,j)
5100     continue
5000   continue
c       write(*,*)'matr=',((matr(i,j),i=1,10),j=1,10)
       call lucmp(kdtr+1,matr,index,vt1) 
       call lusub(kdtr+1,matr,dx,index)
c       stop      
cccccc Use the iterative scheme to improve the accuracey and make ready for
cccccc the later fast calculation, convert the result to real space 
cccccc coordinate. ccccccccccccccccccccccccccccccccccccccccc   
       call torusimpr(1,ktrunc,d,kdtr,dx,dx0,matr,matr0,index,vt1,ct2)
       m=2*trunc+2
       do 7000 i=1,d,1
         do 7300 j=1,2*trunc,1
           datk1(j)=0
7300     continue
         datk1(1)=a(i)+dh*dx(i)
         datk1(2)=a(i+d)+dh*dx(i+d)
         do 7100 j=3,ktrunc,1
           m1=(j-1)*d+i
           t1=a(m1)+dh*dx(m1)
           datk1(j)=t1
         if (mod(j,2).eq.1) then
             datk1(m-j)=t1           
           else
             datk1(m-j+2)=-t1
           end if
7100     continue
         isign=1
         call four1(trunc,isign,datk1) 
         if (mod(i,2).eq.-10) then       
           do 7200 j=1,trunc,1
             y((j-1)*d+i)=datk1(2*(j-1)+1)+(j-1)*t0
7200       continue  
         else
           do 7400 j=1,trunc,1
             y((j-1)*d+i)=datk1(2*(j-1)+1)
7400       continue  
         end if
c         write(*,*) 'datk1=',datk1
c         stop
7000   continue
       omg=omg+dh*dx(kdtr+1)
       end
