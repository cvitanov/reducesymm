function aa=MF(a)
%
% 
ac=a(1:2:size(a))+1i*a(2:2:size(a));
ainv=ac;
phi1=angle(ac(1));
r1=abs(ac(1));
for k=2:2:size(ac), % k even case
    phik=angle(ac(k));
    rk=abs(ac(k));
    ainv(k)=rk*abs(cos(phik-k*phi1)) + 1i*(rk*sin(phik-k*phi1));
end
for k=1:2:size(ac), % k odd case
    phik=angle(ac(k));
    rk=abs(ac(k));
    ainv(k)=rk*cos(phik-k*phi1) + 1i*rk*abs(sin(phik-k*phi1));
end
aa=zeros(size(a));
aa(1:2:size(a))=real(ainv);
aa(2:2:size(a))=imag(ainv);

end