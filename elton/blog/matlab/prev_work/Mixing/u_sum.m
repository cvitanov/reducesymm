% u_sum.m

function u = u_sum(tt,rr)
global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N uspec1 uspec2 uspec3 

x = rr(1);
y = rr(2);
z = rr(3);

zvec = [1/2];
zvec2 = [0];
zvec3 = [0];
for mz = 1:Mz-1
    zvec = [zvec; exp(2*pi*i*mz*z/Lz)];
    zvec2 = [zvec2;(2*pi*i*mz/Lz)*exp(2*pi*i*mz*z/Lz)];
    zvec3 = [zvec3;(-4*(pi*mz/Lz)^2)*exp(2*pi*i*mz*z/Lz)];
end

xvec = [];
xvec2 = [];
xvec3 = [];
for mx = 0:Mx-1
    if 0 <= mx & mx <= Mx/2
        kx = mx;
    else
        kx = mx - Mx;
    end      
    xvec = [xvec; exp(2*pi*i*kx*x/Lx)];
    xvec2 = [xvec2;(2*pi*i*kx/Lx)*exp(2*pi*i*kx*x/Lx)];
    xvec3 = [xvec3;(-4*(pi*kx/Lx)^2)*exp(2*pi*i*kx*x/Lx)];
end

T = zeros(My,1);                    % should generalize to [a,b].
T(1) = 1;        % 0th cheby poly of 1st kind
T(2) = y;
U = zeros(My,1);                    % should generalize to [a,b].
U(1) = 1;        % 0th cheby poly of 2nd kind
U(2) = 2*y;
derT = zeros(My,1);  % derivative of T
derT(1) = 0;
derT(2) = 1;
der2T = zeros(My,1);  % second derivative of T
der2T(1) = 0;
der2T(3) = 0;
for my = 3:My
    T(my) = 2*y*T(my-1) - T(my-2);
    U(my) = 2*y*U(my-1) - U(my-2);
    derT(my) = (my-1)*U(my-1);           % note index shift necessary
%     if y == 1
%         der2T(my) = (1/3)*((my-1)^4 - (my-1)^2);
%     elseif y == -1
%         der2T(my) = ((-1)^(my-1))*(1/3)*((my-1)^4 - (my-1)^2);
%     else    
%         der2T(my) = (my-1)*((my-1)*T(my)-y*U(my-1))/(y^2-1);
%     end
end

v = kron(xvec,zvec);
w = kron(T,v);

u1 = 2*real(uspec1.'*w)+y;       % u1+y adds the laminar flow to u
u2 = 2*real(uspec2.'*w);
u3 = 2*real(uspec3.'*w);
utemp = [u1; u2; u3];  

 v1a = kron(xvec2,zvec);
 v1b = kron(T,v1a);
 v2b = kron(derT,v);
 v3a = kron(xvec,zvec2);
 v3b = kron(T,v3a);
% 
% vv1a = kron(xvec3,zvec);
% vv1b = kron(T,vv1a);
% vv2b = kron(der2T,v);
% vv3a = kron(xvec,zvec3);
% vv3b = kron(T,vv3a);
% 
 A11 = 2*real(uspec1.'*v1b);
 A12 = 2*real(uspec1.'*v2b)+1 ;   % laminar addition
 A13 = 2*real(uspec1.'*v3b);
 A21 = 2*real(uspec2.'*v1b);
 A22 = 2*real(uspec2.'*v2b);
 A23 = 2*real(uspec2.'*v3b); 
 A31 = 2*real(uspec3.'*v1b);
 A32 = 2*real(uspec3.'*v2b);
 A33 = 2*real(uspec3.'*v3b);
% 
% L11 = 2*real(uspec1.'*vv1b);
% L12 = 2*real(uspec1.'*vv2b);
% L13 = 2*real(uspec1.'*vv3b);
% L21 = 2*real(uspec2.'*vv1b);
% L22 = 2*real(uspec2.'*vv2b);
% L23 = 2*real(uspec2.'*vv3b); 
% L31 = 2*real(uspec3.'*vv1b);
% L32 = 2*real(uspec3.'*vv2b);
% L33 = 2*real(uspec3.'*vv3b);
% 
% 
% % matrix of velocity gradients
A = [A11 A12 A13;A21 A22 A23; A31 A32 A33];
%u = [utemp A];
u = utemp;
% trace = A11 + A22 + A33;
% [evecs evals] = eig(A);
% 
% % matrix of Laplacians
% L = [L11 L12 L13;L21 L22 L23; L31 L32 L33];




