function aa=mfinvSO2(a)
%
% 
% 
% dmp=0.1;
% dmp2=0.1;
% ac=a(1:2:size(a))+1i*a(2:2:size(a));
% ainv=ac;
% phi1=angle(ac(1));
% r1=abs(ac(1));
% for k=2:2:size(ac), % k even case
%     phik=angle(ac(k));
%     rk=abs(ac(k));
%     ainv(k)=(1-exp(-r1^2/dmp))*rk*cos(phik-k*phi1) + 1i*(exp(-r1^2/dmp2)*rk^2+(1-exp(-r1^2/dmp))*rk*sin(phik-k*phi1));
% end
% for k=1:2:size(ac), % k odd case
%     phik=angle(ac(k));
%     rk=abs(ac(k));
%     ainv(k)=exp(-r1^2/dmp2)*rk^2+(1-exp(-r1^2/dmp))*rk*cos(phik-k*phi1) + 1i*(1-exp(-r1^2/dmp))*rk*sin(phik-k*phi1);
% end
% aa=zeros(size(a));
% aa(1:2:size(aa))=real(ainv);
% aa(2:2:size(aa))=imag(ainv);

M=4;
ac=a(1:2:size(a))+1i*a(2:2:size(a));
aa=zeros(M*size(a,1),1);
dmp=0.2;
for m=1:M,
    ainv=zeros(size(ac));
    phim=angle(ac(m))/m;
    rm=abs(ac(m));
    for k=m:m:size(ac), 
        phik=angle(ac(k));
        rk=abs(ac(k));
        ainv(k)=(1-exp(-rm^2/dmp))*rk*cos(phik-k*phim) + 1i*((1-exp(-rm^2/dmp))*rk*sin(phik-k*phim));
    end
    aa(2*(m-1)+1:2*M:size(aa))=real(ainv);
%     aa(2*(m-1)+1)=0;
    aa(2*m:2*M:size(aa))=imag(ainv);
end


% m=3;
% dmp=0.1;
% dmp2=0.1;
% ac=a(1:2:size(a))+1i*a(2:2:size(a));
% ainv=ac;
% phim=angle(ac(m))/m;
% rm=abs(ac(m));
% for k=1:size(ac), 
%     phik=angle(ac(k));
%     rk=abs(ac(k));
%     ainv(k)=(1-exp(-rm^2/dmp))*rk*cos(phik-k*phim) + 1i*((1-exp(-rm^2/dmp))*rk*sin(phik-k*phim));
% end
% aa=zeros(size(a));
% aa(1:2:size(aa))=real(ainv);
% aa(2:2:size(aa))=imag(ainv);

% dmp=0.1;
% dmp2=0.1;
% ac=a(1:2:size(a))+1i*a(2:2:size(a));
% ainv=ac;
% r=zeros(size(ac));
% phi1=angle(ac(1));
% r(1)=abs(ac(1));
% for k=2:2:size(ac), % k even case
%     phik=angle(ac(k));
%     r(k)=abs(ac(k));
%     ainv(k)=(1-exp(-r(1)^2/dmp))*r(k)*cos(phik-k*phi1) + 1i*((1-exp(-r(1)^2/dmp))*r(k)*sin(phik-k*phi1));
% end
% for k=1:2:size(ac), % k odd case
%     phik=angle(ac(k));
%     r(k)=abs(ac(k));
%     ainv(k)=(1-exp(-r(1)^2/dmp))*r(k)*cos(phik-k*phi1) + 1i*(1-exp(-r(1)^2/dmp))*r(k)*sin(phik-k*phi1);
% end
% aa=zeros(floor(3*size(a)/2));
% aa(1:3:size(aa))=real(ainv);
% aa(2:3:size(aa))=imag(ainv);
% aa(3:3:size(aa))=exp(-r(1)^2/dmp2)*r.^2;
% aa(3)=0;

% 
% dmp=0.1;
% dmp2=0.1;
% ac=a(1:2:size(a))+1i*a(2:2:size(a));
% ainv=ac;
% phi1=angle(ac(1));
% r1=abs(ac(1));
% for k=2:2:size(ac), % k even case
%     phik=angle(ac(k));
%     rk=abs(ac(k));
%     ainv(k)=(1-exp(-r1^2/dmp))*rk*cos(phik-k*phi1) + 1i*((1-exp(-r1^2/dmp))*rk*sin(phik-k*phi1));
% end
% for k=1:2:size(ac), % k odd case
%     phik=angle(ac(k));
%     rk=abs(ac(k));
%     ainv(k)=(1-exp(-r1^2/dmp))*rk*cos(phik-k*phi1) + 1i*(1-exp(-r1^2/dmp))*rk*sin(phik-k*phi1);
% end
% aa=zeros(size(a));
% aa(1:2:size(aa))=real(ainv);
% aa(2:2:size(aa))=imag(ainv);

% dmp=0.1;
% dmp2=0.1;
% ac=a(1:2:size(a))+1i*a(2:2:size(a));
% ainv=ac;
% phi1=angle(ac(1));
% r1=abs(ac(1));
% for k=2:2:size(ac), % k even case
%     phik=angle(ac(k));
%     rk=abs(ac(k));
%     ainv(k)=(1-exp(-r1^2/dmp))*rk*cos(phik-k*phi1) + 1i*(exp(-r1^2/dmp2)*rk^3+(1-exp(-r1^2/dmp))*rk*sin(phik-k*phi1));
% end
% for k=1:2:size(ac), % k odd case
%     phik=angle(ac(k));
%     rk=abs(ac(k));
%     ainv(k)=exp(-r1^2/dmp2)*rk^3+(1-exp(-r1^2/dmp))*rk*cos(phik-k*phi1) + 1i*(1-exp(-r1^2/dmp))*rk*sin(phik-k*phi1);
% end
% aa=zeros(size(a));
% aa(1:2:size(aa))=real(ainv);
% aa(2:2:size(aa))=imag(ainv);




% dmp=0.1;
% dmp2=1.;
% ac=a(1:2:size(a))+1i*a(2:2:size(a));
% ainv=ac;
% phi1=angle(ac(1));
% r1=abs(ac(1));
% for k=2:2:size(ac), % k even case
%     phik=angle(ac(k));
%     rk=abs(ac(k));
%     ainv(k)=(1-exp(-r1^2/dmp))*rk*cos(phik-k*phi1) + 1i*(exp(-r1^2/dmp2)*rk^2+(1-exp(-r1^2/dmp))*rk*sin(phik-k*phi1));
% end
% for k=1:2:size(ac), % k odd case
%     phik=angle(ac(k));
%     rk=abs(ac(k));
%     ainv(k)=exp(-r1^2/dmp2)*rk^2+(1-exp(-r1^2/dmp))*rk*cos(phik-k*phi1) + 1i*(1-exp(-r1^2/dmp))*rk*sin(phik-k*phi1);
% end
% aa=zeros(size(a));
% aa(1:2:size(aa))=real(ainv);
% aa(2:2:size(aa))=imag(ainv);
end