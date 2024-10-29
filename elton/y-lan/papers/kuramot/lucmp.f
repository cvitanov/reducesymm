c   This routine decomposes matrix a(n,n) into lower and upper triangular matrix
c   which will be stored in a(n,n). Partial pivoting is used so a vector 
c   index(n) has to be used to record the row interchange sequence. It can be 
c   used to solve linear equations and to invert a matrix. 



      subroutine lucmp(n,a,index,vv)
      integer n,index(n)
      double precision a(n,n),vv(n)
      integer i,j,k,imax
      double precision big,dum,sum,temp
            
      do 1000 i=1,n
        big=0
        do 1100 j=1,n
          temp=dabs(a(i,j))
          if(temp.gt.big) big=temp     
1100    continue
        if(big.eq.0) go to 100
        vv(i)=1./big
1000  continue                
      
      do 2000 j=1,n
        do 2100 i=1,j-1
          sum=a(i,j) 
          do 2150 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
2150      continue
          a(i,j)=sum
2100    continue
        big=0.0
        do 2200 i=j,n
          sum=a(i,j)
          do 2250 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
2250      continue
          a(i,j)=sum
          dum=dabs(sum)*vv(i)
          if(dum.ge.big) then
            big=dum
            imax=i
          end if
2200    continue
        if(j.ne.imax) then
          do 2300 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
2300      continue
          vv(imax)=vv(j)
        end if
        index(j)=imax
        if(a(j,j).eq.0) go to 100
        if(j.ne.n) then
          dum=1/a(j,j)
          do 2400 i=j+1,n
            a(i,j)=a(i,j)*dum
2400      continue
        end if
2000  continue                        
      go to 200                                                                         
100   write(*,*) 'Singular matrix in routine lucmp'
200   end