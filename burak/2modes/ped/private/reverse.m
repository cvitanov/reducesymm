function rJ=reverse(J)

n=size(J,1); m=size(J,2)/n; rJ=zeros(n,m*n);
for i=1:m
    rJ(:,(i-1)*n+1:i*n) = J(:,(m-i)*n+1:(m-i+1)*n);
end