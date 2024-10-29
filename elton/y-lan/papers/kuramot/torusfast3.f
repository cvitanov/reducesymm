c This routine is using variational method to find 2-d tori in a general flow 
c provided that the initial guess is quite close to a true torus.

      subroutine torusfast3(omg1,omg2,x,xx,jb,y,ct2,dh,t8,
     \flag3,matr0,index,count3)
      implicit none
      integer d,trunc,dtr,m1,m2,m3,m4,m5,m22,m33,ct2,ktrunc,kdtr,nn
      integer ktr,ii,m0,mm,m6,m7,m8,m9,m11,m44,m55,tr,flag,flag3
      parameter(d=4,trunc=32,nn=2,tr=trunc**nn,dtr=d*tr)
      parameter(ktrunc=trunc/2)
      parameter(ktr=(ktrunc-1)**nn,kdtr=d*ktr)
      integer i,j,k,m,n,isign,index(kdtr+nn),indt(2)
      integer k1,k2,l1,l2,isign1,isign2,count3
      double precision omg1,omg2,x(dtr),xx(dtr),jb(d,d,tr),y(dtr)
      double precision jbfr(d,d,ktr),matr(kdtr+nn,kdtr+nn)
      double precision tempk1(trunc),tempk2(trunc)
      double precision t1,t2,t3,t4,t5,t6,t8,dh,pi,t0
      double precision a(kdtr),b(kdtr),dx(kdtr+nn)
      double precision datk1(2*trunc),datk2(2*trunc)
      double precision dx0(kdtr+nn),vt1(kdtr+nn)
      double precision z(dtr),zz(dtr)
      double precision matr0(kdtr+nn,kdtr+nn),aa(kdtr),bb(kdtr)
      double precision temp1((ktrunc-1)*d,trunc)
      double precision temp2((ktrunc-1)*d,trunc)
      double precision temp3((ktrunc-1)*d,ktrunc-1)
      double precision temp4((ktrunc-1)*d,ktrunc-1)
     
      pi=3.1415926535897932D0
      t0=2*pi/trunc
      n=ktrunc/2-1
      do 500 i=1,kdtr+nn,1
        dx(i)=0
        do 510 j=1,kdtr+nn,1
          matr(i,j)=0
510     continue
c        if(mod(i,2).eq.0.and.i.le.10) matr(kdtr+1,i)=1
500   continue

cccccc From x,xx,jb, generate the Fourier coefficients a,b,jbfr ccccccccccccc
      flag=1
      call fint2(d,trunc,tr,ktrunc,dtr,kdtr,tempk1,tempk2,
     \datk1,datk2,temp1,temp2,temp3,temp4,x,xx,a,b,t0,flag)     
      do 2000 i=1,d/2,1
        m1=2*(i-1) 
        do 2100 j=1,tr,1
          m2=(j-1)*d
          do 2110 k=1,d,1
            m3=m2+k
            z(m3)=jb(k,m1+1,j)
            zz(m3)=jb(k,m1+2,j)
2110      continue
2100    continue
        flag=0
        call fint2(d,trunc,tr,ktrunc,dtr,kdtr,tempk1,tempk2,
     \datk1,datk2,temp1,temp2,temp3,temp4,z,zz,aa,bb,t0,flag)
        do 2200 j=1,ktr,1
          m2=(j-1)*d
          do 2210 k=1,d,1
            m3=m2+k
            jbfr(k,m1+1,j)=aa(m3)
            jbfr(k,m1+2,j)=bb(m3)
2210      continue
2200    continue
2000  continue
      write(*,*) 'b value=',(b(i),i=1,16)
      write(*,*) 'xx value=',(xx(i),i=1,16)
c      stop
cccccc Construct the big matrix matr in the variational equation ccccccccccc
       write(*,*) 'construct'
       do 4000 k1=0,n,1
         m0=k1*(ktrunc-1)*d*2
         if(k1.ne.0) m0=m0-(ktrunc-1)*d
         m1=m0+(ktrunc-1)*d
         do 4100 k2=0,n,1
           m2=m0+k2*d*2           
           if(k2.ne.0) m2=m2-d
           m3=m2+d
           t1=cos(k1*omg1+k2*omg2)
           t2=sin(k1*omg1+k2*omg2)
           do 4110 i=1,d,1  
             m22=m2+i
             if (k2.eq.0) then
               m33=m1+i
             else
               m33=m3+i
             end if 
             if (k2.eq.0.and.k1.eq.0) then
               matr(m22,m22)=t1
               dx(m22)=b(m22)-a(m22)*t1
               if(i.eq.1) then
                 matr(m22,kdtr+1)=matr(m22,kdtr+1)+1
                 dx(m22)=dx(m22)-omg1
               end if
               if(i.eq.3) then
                 matr(m22,kdtr+2)=matr(m22,kdtr+2)+1
                 dx(m22)=dx(m22)-omg2              
               end if
             else  
               matr(m22,m22)=t1
               matr(m22,m33)=-t2
               matr(m33,m22)=t2
               matr(m33,m33)=t1
               matr(m22,kdtr+1)=-k1*(a(m22)*t2+a(m33)*t1)
               matr(m33,kdtr+1)=k1*(a(m22)*t1-a(m33)*t2) 
               matr(m22,kdtr+2)=-k2*(a(m22)*t2+a(m33)*t1)
               matr(m33,kdtr+2)=k2*(a(m22)*t1-a(m33)*t2)
               dx(m22)=b(m22)-a(m22)*t1+a(m33)*t2
               dx(m33)=b(m33)-a(m22)*t2-a(m33)*t1 
             end if
4110       continue
           do 3000 l1=k1-n,n,1
c             m4=abs(l1)*(ktrunc-1)*d*2
c             if(l1.ne.0) m4=m4-(ktrunc-1)*d
             do 3100 l2=k2-n,n,1 
               m4=abs(l1)*(ktrunc-1)*d*2
               if(l1.ne.0) m4=m4-(ktrunc-1)*d
               isign=1
               isign1=1
               isign2=1
c               m5=m4+abs(l2)*d*2
               if (l2.lt.0) then
                 isign2=-1
               elseif(l2.eq.0.and.l1.lt.0) then
                 isign2=-1
               end if
               m6=k2-l2
               if (m6.lt.0) then
                 isign1=-1
               elseif(m6.eq.0.and.k1-l1.lt.0) then
                 isign1=-1
               end if
               m7=l1*isign2 
               m8=(k1-l1)*isign1
               if(m7.lt.0.and.l2.ne.0) m4=m4+(ktrunc-1)*d
               m4=m4+abs(l2)*d*2
               if(l2.ne.0) m4=m4-d
               m9=abs(m8)*(ktrunc-1)*2
               if(m8.ne.0) m9=m9-(ktrunc-1)
               if(m8.lt.0.and.m6.ne.0) m9=m9+ktrunc-1  
               m9=m9+abs(m6)*2+1
               if(m6.ne.0) m9=m9-1
               if (m6.eq.0) then
                 m5=m9+ktrunc-1
               else
                 m5=m9+1
               end if
               do 3110 i=1,d,1
                 m22=m2+i
                 if (k2.eq.0) then
                   m33=m1+i
                 else
                   m33=m3+i
                 end if 
                 do 3111 j=1,d,1
                   m44=m4+j
                   if (l2.eq.0) then
                     m55=m44+(ktrunc-1)*d
                   else 
                     m55=m44+d
                   end if
         matr(m22,m44)=matr(m22,m44)-jbfr(i,j,m9)
         if(m4.ne.0.and.m9.ne.1) matr(m22,m55)=matr(m22,m55)
     \+isign1*isign2*jbfr(i,j,m5)
         if (m2.ne.0) then
         if(m9.ne.1) matr(m33,m44)=matr(m33,m44)-isign1*jbfr(i,j,m5)
         if(m4.ne.0) matr(m33,m55)=matr(m33,m55)-isign2*jbfr(i,j,m9)      
         end if
c         if (matr(13,117)+1.lt.1E-3) then
c           write(*,*) 'jbfr=',jbfr(i,j,m5),'m22,m33,m44,m55,m4,m5,m6,
c     \m7,m8,m9,l1,l2,k1,k2'
c           write(*,*) m22,m33,m44,m55,m4,m5,m6,m7,m8,m9,l1,l2,k1,k2
c           stop 
c         end if
3111             continue
3110           continue
3100         continue
3000       continue
4100     continue
c      write(*,*) 'matr value=',(matr(13,j),j=115,120)
c      stop
         if (k1.ne.0) then
         do 4200 k2=1,n,1
           m2=m1+k2*d*2-d           
           m3=m2+d
           t1=cos(-k1*omg1+k2*omg2)
           t2=sin(-k1*omg1+k2*omg2)
           do 4210 i=1,d,1  
             m22=m2+i
             m33=m3+i
             matr(m22,m22)=t1
             matr(m22,m33)=-t2
             matr(m33,m22)=t2
             matr(m33,m33)=t1
             matr(m22,kdtr+1)=k1*(a(m22)*t2+a(m33)*t1)
             matr(m33,kdtr+1)=-k1*(a(m22)*t1-a(m33)*t2) 
             matr(m22,kdtr+2)=-k2*(a(m22)*t2+a(m33)*t1)
             matr(m33,kdtr+2)=k2*(a(m22)*t1-a(m33)*t2)
             dx(m22)=b(m22)-a(m22)*t1+a(m33)*t2
             dx(m33)=b(m33)-a(m22)*t2-a(m33)*t1 
4210       continue
           do 3500 l1=-n,-k1+n,1
c             m4=abs(l1)*(ktrunc-1)*d*2
c             if(l1.ne.0) m4=m4-(ktrunc-1)*d
             do 3510 l2=k2-n,n,1 
               m4=abs(l1)*(ktrunc-1)*d*2
               if(l1.ne.0) m4=m4-(ktrunc-1)*d
               isign=1
               isign1=1
               isign2=1
c               m5=m4+abs(l2)*d*2
               if (l2.lt.0) then
                 isign2=-1
               elseif(l2.eq.0.and.l1.lt.0) then
                 isign2=-1
               end if
               m6=k2-l2
               if (m6.lt.0) then
                 isign1=-1
               elseif(m6.eq.0.and.-k1-l1.lt.0) then
                 isign1=-1
               end if
               m7=l1*isign2 
               m8=(-k1-l1)*isign1
               if(m7.lt.0.and.l2.ne.0) m4=m4+(ktrunc-1)*d
               m4=m4+abs(l2)*d*2
               if(l2.ne.0) m4=m4-d
               m9=abs(m8)*(ktrunc-1)*2
               if(m8.ne.0) m9=m9-(ktrunc-1)
               if(m8.lt.0.and.m6.ne.0) m9=m9+ktrunc-1  
               m9=m9+abs(m6)*2+1
               if(m6.ne.0) m9=m9-1
               if (m6.eq.0) then
                 m5=m9+ktrunc-1
               else
                 m5=m9+1
               end if
               do 3511 i=1,d,1
                 m22=m2+i
                 m33=m3+i 
                 do 3512 j=1,d,1
                   m44=m4+j
                   if (l2.eq.0) then
                     m55=m44+(ktrunc-1)*d
                   else 
                     m55=m44+d
                   end if
         matr(m22,m44)=matr(m22,m44)-jbfr(i,j,m9)
         if(m4.ne.0.and.m9.ne.1) matr(m22,m55)=matr(m22,m55)
     \+isign1*isign2*jbfr(i,j,m5)
         if(m9.ne.1) matr(m33,m44)=matr(m33,m44)-isign1*jbfr(i,j,m5)
         if(m4.ne.0) matr(m33,m55)=matr(m33,m55)-isign2*jbfr(i,j,m9)
3512             continue
3511           continue
3510         continue
3500       continue
4200     continue
         end if
4000   continue
c       do 4500 i=1,kdtr+nn,1
c         write(*,*)'i=',i
c         write(*,*)'ccccccccccccccccccccccccc'
c         do 4510 j=1,kdtr+nn,1
c           write(*,'(F12.6)') matr(i,j)
c4510     continue
c4500   continue
c       write(*,*) 't8=',t8, ' ind=',indt
c       stop 
       matr(kdtr+1,kdtr+1)=1
       matr(kdtr+2,kdtr+2)=1

cccccc Solve the linear equation for the change of the Fourier coeffient cccc
       write(*,*) 'solve'
       t8=0
       if(flag3.eq.0) then 
       flag3=1
       count3=1
c       t8=0
       mm=0
       do 5000 i=1,kdtr+nn,1
         t1=dx(i)
         dx0(i)=t1
         if(t8.lt.abs(t1)) mm=i
         t8=max(t8,abs(t1))
         do 5100 j=1,kdtr+nn,1
           matr0(i,j)=matr(i,j)
5100     continue
5000   continue
       call lucmp(kdtr+nn,matr,index,vt1) 
       call lusub(kdtr+nn,matr,dx,index)
c       write(*,*) 'The biggest difference mm=',mm,t8
c       write(*,*) 'jbfr value=',(((jbfr(i,j,k),i=1,d),j=1,d),k=1,2)
c       write(*,*) 'matr0 value=',(matr0(13,j),j=115,120)
c       stop
      
cccccc Use the iterative scheme to improve the accuracey and make ready for
cccccc the later fast calculation, convert the result to real space 
cccccc coordinate. ccccccccccccccccccccccccccccccccccccccccc  
       write(*,*) 'improve' 
       call torusimpr3(nn,kdtr,dx,dx0,matr,matr0,
     \index,vt1,ct2)
       do 5500 i=1,kdtr+nn,1
         do 5510 j=1,kdtr+nn,1
           matr0(i,j)=matr(i,j)
5510     continue
5500   continue
       else
       do 5600 i=1,kdtr+nn,1
         t1=dx(i)
         dx0(i)=t1
         if(t8.lt.abs(t1)) mm=i
         t8=max(t8,abs(t1))
5600   continue
       count3=count3+1
       call lusub(kdtr+nn,matr0,dx,index)
       call torusimpr3(nn,kdtr,dx,dx0,matr0,matr,
     \index,vt1,ct2)
       write(*,*) 'The biggest difference mm=',mm,t8
       end if 
       m=2*trunc+1
c       write(*,*) 'dx=',(dx(i),i=1,16,1)
       t1=0
       mm=0
       do 6000 i=1,kdtr,1
         a(i)=a(i)+dh*dx(i)
c         if(t1.lt.abs(dx(i))) mm=i         
c         t1=max(t1,abs(dx(i)))
6000   continue
c       write(*,*) 'max dx',t1,i,dh
       do 4600 i=d+1,(ktrunc-1)*d,1
         temp3(i,1)=a(i)
4600   continue
       do 4500 i=1,ktrunc-1,1
         do 4510 j=1,d,1
           m1=j+(i-1)*(ktrunc-1)*d
           temp3(j,i)=a(m1)
4510     continue
4500   continue
       do 7000 i=1,ktrunc/2-1,1
         m1=(i-1)*(ktrunc-1)*d*2+(ktrunc-1)*d
         m11=m1+(ktrunc-1)*d 
         m3=i*2
         m33=m3+1
         do 7100 j=1,ktrunc/2-1,1
           m2=(j-1)*d*2+d
           m22=m2+d            
           do 7110 k=1,d,1
             temp3(m2+k,m3)=(a(m1+m2+k)+a(m11+m2+k))/2.
             temp3(m2+k,m33)=(a(m1+m22+k)-a(m11+m22+k))/2.
             temp3(m22+k,m3)=(a(m1+m22+k)+a(m11+m22+k))/2.
             temp3(m22+k,m33)=(a(m11+m2+k)-a(m1+m2+k))/2.            
7110       continue
7100     continue
7000   continue
       do 8000 i=1,(ktrunc-1)*d,1
         do 8100 j=1,2*trunc,1
           datk1(j)=0
8100     continue
         datk1(1)=temp3(i,1)
         datk1(2)=0
         do 8200 j=2,ktrunc-1,1
           t1=temp3(i,j)
           datk1(j+1)=t1
         if (mod(j,2).eq.0) then
             datk1(m-j)=t1           
           else
             datk1(m-j+2)=-t1
           end if
8200     continue
         isign=1
         call four1(trunc,isign,datk1) 
         do 8300 j=1,trunc,1
           temp1(i,j)=datk1(2*(j-1)+1)
8300     continue
8000   continue
       do 8500 i=1,trunc,1
       m2=(i-1)*(ktrunc-1)*d
       m22=(i-1)*trunc*d
       do 8510 k=1,d,1
         do 8511 j=1,2*trunc,1
           datk1(j)=0
8511     continue
         datk1(1)=temp1(k,i)
         datk1(2)=0
         do 8512 j=2,ktrunc-1,1
           m1=(j-1)*d+k
           t1=temp1(m1,i)
           datk1(j+1)=t1
         if (mod(j,2).eq.0) then
             datk1(m-j)=t1           
           else
             datk1(m-j+2)=-t1
           end if
8512     continue
         isign=1
         call four1(trunc,isign,datk1)         
         do 8513 j=1,trunc,1
           y(m22+(j-1)*d+k)=datk1(2*(j-1)+1)
8513     continue
8510   continue
8500   continue
       do 9000 i=1,trunc,1
          m1=(i-1)*d*trunc
          t1=(i-1)*t0
          do 9100 j=1,trunc,1
            m2=m1+(j-1)*d
            t2=(j-1)*t0
            y(m2+1)=y(m2+1)+t1
            y(m2+3)=y(m2+3)+t2
9100      continue
9000   continue
       omg1=omg1+dh*dx(kdtr+1)
       omg2=omg2+dh*dx(kdtr+2)
       write(*,*) 'back'
       end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCThe following is used to calculate two dimensional Fourier modesCCCCC
CCCCcoefficients CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fint2(d,trunc,tr,ktrunc,dtr,kdtr,tempk1,tempk2,
     \datk1,datk2,temp1,temp2,temp3,temp4,x,xx,a,b,t0,flag)
      integer i,j,k,d,trunc,tr,ktrunc,dtr,kdtr
      integer m1,m2,m3,flag,m11,m22,m1r,m1i
      double precision tempk1(trunc),tempk2(trunc)
      double precision datk1(2*trunc),datk2(2*trunc)
      double precision temp1((ktrunc-1)*d,trunc)
      double precision temp2((ktrunc-1)*d,trunc)
      double precision x(dtr),xx(dtr),a(kdtr),b(kdtr) 
      double precision t0,t1,t2,temp3((ktrunc-1)*d,ktrunc-1)
      double precision temp4((ktrunc-1)*d,ktrunc-1)

      if(flag.eq.1) then
        do 500 i=1,trunc,1
          m1=(i-1)*d*trunc
          t1=(i-1)*t0
          do 510 j=1,trunc,1
            m2=m1+(j-1)*d
            t2=(j-1)*t0
            x(m2+1)=x(m2+1)-t1
            x(m2+3)=x(m2+3)-t2
            xx(m2+1)=xx(m2+1)-t1
            xx(m2+3)=xx(m2+3)-t2
510       continue
500     continue
      end if
      do 900 k=1,trunc,1
      m2=(k-1)*d*trunc
      do 1000 i=1,d,1
        do 1100 j=1,trunc,1
          m1=m2+(j-1)*d+i
c          t1=(j-1)*t0
c          if (i.eq.3.and.flag.eq.1) then
c            tempk1(j)=x(m1)-t1
c            tempk2(j)=xx(m1)-t1
c          else
            tempk1(j)=x(m1)
            tempk2(j)=xx(m1)
c          end if
1100    continue
        call dftint(trunc,tempk1,tempk2,datk1,datk2)
        temp1(i,k)=tempk1(1)
        temp2(i,k)=tempk2(1)
        do 1200 j=3,ktrunc,1
          m1=(j-2)*d+i
          temp1(m1,k)=tempk1(j)
          temp2(m1,k)=tempk2(j)
1200    continue 
1000  continue
900   continue
      do 1500 k=1,ktrunc-1,1
      m2=(k-1)*d
      do 1510 i=1,d,1
        m1=m2+i
        do 1511 j=1,trunc,1
c          t1=(j-1)*t0
c          if (i.eq.1.and.flag.eq.1) then
c            tempk1(j)=temp1(m1,j)-t1
c            tempk2(j)=temp2(m1,j)-t1
c          else
            tempk1(j)=temp1(m1,j)
            tempk2(j)=temp2(m1,j)
c          end if
1511    continue
        call dftint(trunc,tempk1,tempk2,datk1,datk2)
        temp3(m1,1)=tempk1(1)
        temp4(m1,1)=tempk2(1)
        do 1512 j=3,ktrunc,1
          temp3(m1,j-1)=tempk1(j)
          temp4(m1,j-1)=tempk2(j)
1512    continue
1510  continue
1500  continue
      do 2000 i=1,ktrunc/2-1,1
        m1=(i-1)*(ktrunc-1)*d*2+(ktrunc-1)*d
        m11=m1+(ktrunc-1)*d
        m1r=i*2
        m1i=m1r+1
        do 2100 j=1,ktrunc/2-1,1
          m2=(j-1)*d*2+d
          m22=m2+d
          do 2110 k=1,d,1
            a(m1+m2+k)=temp3(m2+k,m1r)-temp3(m22+k,m1i)
            a(m1+m22+k)=temp3(m2+k,m1i)+temp3(m22+k,m1r)
            a(m11+m2+k)=temp3(m2+k,m1r)+temp3(m22+k,m1i)
            a(m11+m22+k)=temp3(m22+k,m1r)-temp3(m2+k,m1i)
            b(m1+m2+k)=temp4(m2+k,m1r)-temp4(m22+k,m1i)
            b(m1+m22+k)=temp4(m2+k,m1i)+temp4(m22+k,m1r)
            b(m11+m2+k)=temp4(m2+k,m1r)+temp4(m22+k,m1i)
            b(m11+m22+k)=temp4(m22+k,m1r)-temp4(m2+k,m1i)
2110      continue
2100    continue
2000  continue
      do 3000 i=d+1,(ktrunc-1)*d,1
        a(i)=temp3(i,1)
        b(i)=temp4(i,1) 
3000  continue
      do 4000 i=1,ktrunc-1,1
        do 4100 j=1,d,1
          m1=j+(i-1)*(ktrunc-1)*d
          a(m1)=temp3(j,i)
          b(m1)=temp4(j,i)
4100    continue
4000  continue
      if(flag.eq.1) then
        do 5000 i=1,trunc,1
          m1=(i-1)*d*trunc
          t1=(i-1)*t0
          do 5100 j=1,trunc,1
            m2=m1+(j-1)*d
            t2=(j-1)*t0
            x(m2+1)=x(m2+1)+t1
            x(m2+3)=x(m2+3)+t2
            xx(m2+1)=xx(m2+1)+t1
            xx(m2+3)=xx(m2+3)+t2
5100      continue
5000    continue
      endif
      end
