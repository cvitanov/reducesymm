function [rm, s] = manif2rm(ao, a, nhit, sgn)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
 
if nargin<2, sgn=1; end;
 
[a1, ds1] = keepCont(a(:,1:nhit:end),0.1); % Eliminate need to provide cutoff
[a2, ds2] = keepCont(a(:,2:nhit:end),0.1);
 
s=ds2s([norm(ao-a1(:,1)) ds1 norm(a1(:,end)-a2(:,1)) ds2, sgn]);

rm= s2rm(s(2:end), [size(a1,2), size(a2,2)]);
 
rm=[[0 0]; rm];
 
 
end

