function [cn,xn]=normc(x);
% [cn,xn]=normc(x); gives a row vector of column 2 norms of x
[n p]=size(x);
if issparse(x), xn=sparse(n,p); else xn=zeros(n,p); end;
cn=zeros(1,p);
for j=1:p, cn(j)=norm(x(:,j)); xn(:,j)=x(:,j)/cn(j); end;
