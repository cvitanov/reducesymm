function [ s ] = ds2s(ds, sgn)
% Sum distance of points to parametrize a curve by Euclidian distance.
% First element corresponds to 0. Allow for negative s (set optional
% argument sgn to -1).  

if nargin<2, sgn=1; end;

s=zeros(length(ds)+1,1);
for j=2:size(s),
    s(j)=s(j-1)+sgn*ds(j-1);
end

