function sout  = p2s(a, ao, aman, s, cut )
% Map point to s on an unstable manifold
%   Detailed explanation goes here

aman1=[ao,aman];

d=zeros(size(aman1,2),1);

for i=1:size(aman1,2)
    d(i)=norm(aman1(:,i)-a);
end

dmin=min(d);

if dmin < cut
    posmin = find(d==dmin);
    sout=s(posmin);
else
    sout=NaN;
end

end

