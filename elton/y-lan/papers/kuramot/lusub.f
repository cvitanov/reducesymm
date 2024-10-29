c   This routine follows lucmp to give the solution of linear equations
c   a*x=b and the result will put in b(n).




      subroutine lusub(n,a,b,index)
      integer n,index(n)
      double precision a(n,n),b(n)
      integer i,j,ii,ip
      double precision sum
      ii=0
      
      do 1000 i=1,n
        ip=index(i)
        sum=b(ip)
        b(ip)=b(i)
        if(ii.ne.0) then
          do 1100 j=ii,i-1
            sum=sum-a(i,j)*b(j)
1100      continue           
        else if(sum.ne.0) then
          ii=i
        end if
        b(i)=sum
1000  continue
      do 2000 i=n,1,-1
        sum=b(i)
        do 2100 j=i+1,n            
          sum=sum-a(i,j)*b(j)
2100    continue
        b(i)=sum/a(i,i)
2000  continue
      end                  