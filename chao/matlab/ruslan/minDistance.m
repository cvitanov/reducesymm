function [dmin, aa, bb]=minDistance(a,b,slc)

 dim=size(a,1)/2;
 aa=a; bb=b;
 aa(:,1)=mf(a(:,1),slc);
 bb(:,1)=mf(b(:,1),slc);

 dmin=norm(aa(:,1)-bb(:,1));
 
 disp(dmin);
 
 
 
 for i=1:size(a,2),
     aa(:,i)=mf(a(:,i),slc);
     for j=1:size(b,2),
         slc1=slc;
         bb(:,j)=mf(b(:,j),slc);
         d=norm(aa(:,i)-bb(:,j));
         while d > (1+1000*eps)*norm(a(:,i)-b(:,j)) && slc1 < dim  ,
             slc1=slc1+1;
             disp(['Suspected jump ', num2str(i), ' ', num2str(j)]); 
             aatmp=mf(a(:,i),slc1);
             bbtmp=mf(b(:,j),slc1);
             d=norm(aatmp-bbtmp);
             if d > (1+1000*eps)*norm(a(:,i)-b(:,j)), 
                 disp(['alternative slice failure, k= ', num2str(slc1)]);
             else
                 disp(['success, k= ', num2str(slc1)]);
             end;
         end;
         dmin=min(d,dmin);
     end
 end
