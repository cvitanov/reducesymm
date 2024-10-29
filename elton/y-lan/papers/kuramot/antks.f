      subroutine antks(L,n0,nu,a,k)
      integer n0,n1,i,j
      double precision a(n0)
      double precision k(n0)
      double precision nu,L
      
      double precision x,q,kq

      q=2*dble(3.1415926536)/L
      n1=n0
              
      do 1000 i = 1,n1,1
        x = 0
        kq=i*q
        k(i) = (kq**2-nu*kq**4)*a(i)
        do 1100 j = 1,n1-i,1
          x=x+a(j)*a(j+i)
1100    continue
        x = 2*x
        do 1200 j=1,i-1,1
          x=x-a(j)*a(i-j)
1200    continue
        k(i)=k(i)+x*(kq)
1000    continue
      end
      
     
      subroutine ksjcb(L,n0,nu,dv,a)
      integer n0,n1,i,j
      double precision dv(n0,n0)
      double precision a(n0)
      double precision nu,L,q,kq

      q=2*dble(3.1415926536)/L
      n1=n0

      do 1000 i=1,n0,1
        do 1100 j=1,n0,1
          dv(i,j)=0
1100    continue
1000  continue

      do 2000 i=1,n1,1
        kq=i*q
        do 2100 j=1,n1-i,1
          dv(i,j)=dv(i,j)+2*kq*a(j+i)
          dv(i,j+i)=dv(i,j+i)+2*kq*a(j)
2100    continue
        do 2200 j=1,i-1,1
          dv(i,j)=dv(i,j)-kq*a(i-j)
          dv(i,i-j)=dv(i,i-j)-kq*a(j)
2200    continue
        dv(i,i)=dv(i,i)+kq**2-nu*kq**4
2000  continue
      end

      subroutine ksjcbj(L,n0,nu,dv,a,jcb)
      integer n0,n1,m,m2,i,j,k1,k
      double precision dv(n0,n0)
      double precision a(n0),jcb(n0,n0)
      double precision nu,L,q,kq,t1

      q=2*dble(3.1415926536)/L
      n1=n0

      do 1000 i=1,n0,1
        do 1100 j=1,n0,1
          dv(i,j)=0
1100    continue
1000  continue

      do 2000 i=1,n1,1
        kq=i*q
        do 2100 j=1,n1-i,1
          dv(i,j)=dv(i,j)+2*kq*a(j+i)
          dv(i,j+i)=dv(i,j+i)+2*kq*a(j)
2100    continue
        do 2200 j=1,i-1,1
          dv(i,j)=dv(i,j)-kq*a(i-j)
          dv(i,i-j)=dv(i,i-j)-kq*a(j)
2200    continue
        dv(i,i)=dv(i,i)+kq**2-nu*kq**4
2000  continue
      do 3000 i=1,n0,1
        do 3100 j=1,n0,1
          a(j)=dv(i,j)
3100    continue
        do 3200 j=1,n0,1
          t1=0
          do 3210 k=1,n0,1
            t1=t1+a(k)*jcb(k,j)
3210      continue
          dv(i,j)=t1
3200    continue
3000  continue
      end
