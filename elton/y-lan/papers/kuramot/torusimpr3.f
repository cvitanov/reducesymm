      subroutine torusimpr3(nn,dtr,dx,dx0,matr,matr0,index,
     \vt1,count)
      integer i,j,itmax,dtr,index(dtr+nn)
      double precision matr(dtr+nn,dtr+nn),dx(dtr+nn),dx0(dtr+nn)
      double precision t1,t2,t3,matr0(dtr+nn,dtr+nn),vt1(dtr+nn) 
      integer count   

      itmax=32
      count=0
      t2=1
      eps=1.e-14
      
      do 500 while(t2.gt.eps.and.count.lt.itmax)
      count=count+1
      t2=0
      do 1000 i=1,dtr+nn,1
        t1=0
        do 1100 j=1,dtr+nn,1
          t1=t1+matr0(i,j)*dx(j)
1100    continue
        t3=dx0(i)-t1 
        vt1(i)=t3
        t2=max(t2,abs(t3))
1000  continue
      call lusub(dtr+nn,matr,vt1,index)
      do 2000 i=1,dtr+nn,1
        dx(i)=dx(i)+vt1(i)
2000  continue
500   continue
      write(*,*) 'The discrepancy=',t2,' itrate=',count
      end
