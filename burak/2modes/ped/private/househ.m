function [a,q]=househ(a,k,flag)
% househ() fulfill onestep household transformation to matrix a at column
%k. q is the unitary transform matrix, a is the transformed matrix.
% if flag==1, a(k+1:n,k)=0; => make the under-diagonal elements be zero.
%otherwise, a(k+2:n,k)=0;  => make the subunder-diagonal elements be zero.

    if nargin<3, flag=1; end  %default value is flag=1.
    m=size(a,1); n=size(a,2);  q=zeros(m,m);
    %if m<n, disp('Household error: dimension m<n'); end
    if k>n || k>m, disp('Household error: k>n or k>m'); return; end
    if flag==1, x=a(k:m,k); e=zeros(m-k+1,1); e(1)=1;
        v=sign(x(1))*norm(x)*e+x; v=v/norm(v);
        q(1:k-1,1:k-1)=eye(k-1); q(k:m,k:m)=eye(m-k+1)-2*v*v'; % form the transform matrix.
        %only the under two matrix blocks need to be transformed.
        a(k:m,k:n)=q(k:m,k:m)*a(k:m,k:n); a(k:m,1:k-1)=q(k:m,k:m)*a(k:m,1:k-1); 
        if k+1<=n, a(k+1:n,k)=0; end  % get rid of rundoff error.
    else
        if k+1>n || k+1>m, disp('Hous_hess error: k+1>n or k+1>m'); return; end
        x=a(k+1:m,k); e=zeros(m-k,1); e(1)=1;
        v=sign(x(1))*norm(x)*e+x; v=v/norm(v);
        q(1:k,1:k)=eye(k); q(k+1:m,k+1:m)=eye(m-k)-2*v*v'; % form the transform matrix
        a(k+1:m,k:n)=q(k+1:m,k+1:m)*a(k+1:m,k:n); a(k+1:m,1:k-1)=q(k+1:m,k+1:m)*a(k+1:m,1:k-1);
        if k+2<=n, a(k+2:n,k)=0; end% get rid of rundoff error.
    end
end