function [ aa ] = ksmf1Refl( a )
% Apply reflection operator to KS point on slice c_1=0. 
aa=a;
for k=1:size(aa)/2,
    aa(2*(k-1)+1) = (-1)^(k+1)*a(2*(k-1)+1);
    aa(2*k) = (-1)^k*a(2*k);
end

