function [T,Y,TE,YE,IE]=hextraj(iC)

%Atomic/Potential parameters
%V   = @(r,Z) Z./(r.^2+1).^0.5;
%dV  = @(r,Z) (r*Z)./(r.^2+1).^(1.5);

rr = sqrt(3)/2;
rh = 1;
vv = 2;
hh = 1;
mu = 0;
vc = vv;
dV = @(R) hh/rr/2*(1-cos(R*2*pi/rr));
dF = @(V) mu*(vc - V); 
V  = @(R) hh/rr/2*(rr-R+rr/2/pi*sin(R*2*pi/rr)).*(R-rr<=0);

%Number of initial conditions
nI = 1;
%Number of dimensions
nD = 4;
%Randomization seeding
sd=100*(sum(clock));
stream=RandStream('mt19937ar','Seed',sd);
RandStream.setGlobalStream(stream);

%Integration tolerances
abstol=1e-9;
reltol=1e-9;
maxstep=1e-2;
options=odeset('AbsTol',abstol,'RelTol',reltol,'MaxStep',maxstep,'Events',@event);

[T,Y,TE,YE,IE]=ode45(@sys,0:1e-2:20,iC,options,nI,nD,dV,dF,rr,rh);

% KE=(Y(:,3).^2+Y(:,4).^2)/2;
% 
% [wx wy]=hexwrap(Y(:,1)', Y(:,2)', rr*2/sqrt(3));
% ty=[Y(:,1)'-wx;Y(:,2)'-wy];
% [th, dr] = cart2pol(ty(1,:),ty(2,:));
% VE=V(dr);
% figure;
% plot(T, KE+VE');

% figure; plot(Y(:,1),Y(:,2),'r',YE(:,1),YE(:,2),'og');


function dy = sys(t,y,nI,nD,dV,dF,rr,rh)
t
dy=zeros(nD*nI,1);
% square mod
% ty=mod(y,2*rr);
% idx=ty>rr;
% ty(idx)=ty(idx)-2*rr;
% triang mod
ty=zeros(nD*nI,1);
[wx wy]=hexwrap(y(1:nD:end)', y(2:nD:end)', rh);
ty(1:nD:end)=y(1:nD:end)-wx';
ty(2:nD:end)=y(2:nD:end)-wy';
%%
% same for all lattice
dr=sqrt(ty(1:nD:end).^2+ty(2:nD:end).^2);
dv=sqrt(y(3:nD:end).^2+y(4:nD:end).^2);
dy(1:nD:end)=y(3:nD:end);
dy(2:nD:end)=y(4:nD:end);
F1=dV(dr);
F2=dF(dv);
dy(3:nD:end)=F1.*ty(1:nD:end)./dr.*(dr-rr<=0)+F2.*y(3:nD:end)./dv;
dy(4:nD:end)=F1.*ty(2:nD:end)./dr.*(dr-rr<=0)+F2.*y(4:nD:end)./dv;


function [value isTerminal direction] = event(~,y,nI,nD,~,~,rr,rh)
% square mod
% ty=mod(y,2*rr);
% idx=ty>rr;
% ty(idx)=ty(idx)-2*rr;
% triang lattice

ty=zeros(nD*nI,1);
[wx wy]=hexwrap(y(1:nD:end)', y(2:nD:end)', rh);
ty(1:nD:end)=y(1:nD:end)-wx';
ty(2:nD:end)=y(2:nD:end)-wy';
% rest is the same
value = sqrt(ty(1:nD:end).^2+ty(2:nD:end).^2) - rr;         %Function of phase space variables that is set equal to 0.   f(x,p) = 0.

direction=0*ones(nI,1);                                     %Direction of crossing
isTerminal=zeros(nI,1);                                     %Does the integrator stop after crossing
