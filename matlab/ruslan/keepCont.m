function [aa, anorm]= keepCont(a, cut)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

icut=size(a,2);

anorm=[];

for i=2:size(a,2),
     if (norm(a(:,i)-a(:,i-1))> cut), 
         icut=i-1; 
         break; 
     end;
     anorm=[anorm,norm(a(:,i)-a(:,i-1))];
end

 aa= a(:,1:icut);    

end

