% u_sum.m

function u = u_sum(tt,rr)
global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N uspec1 uspec2 uspec3 evals

x = rr(1);
y = rr(2);
z = rr(3);

zvec = [1/2];
zvec2 = [0];
for mz = 1:Mz-1
    zvec = [zvec; exp(2*pi*i*mz*z/Lz)];
    zvec2 = [zvec2;(2*pi*i*mz/Lz)*exp(2*pi*i*mz*z/Lz)];
end

xvec = [];
xvec2 = [];
for mx = 0:Mx-1
    if 0 <= mx & mx <= Mx/2
        kx = mx;
    else
        kx = mx - Mx;
    end      
    xvec = [xvec; exp(2*pi*i*kx*x/Lx)];
    xvec2 = [xvec2;(2*pi*i*kx/Lx)*exp(2*pi*i*kx*x/Lx)];
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
for my = 3:My
    T(my) = 2*y*T(my-1) - T(my-2);
    U(my) = 2*y*U(my-1) - U(my-2);
    derT(my) = (my-1)*U(my-1);           % note index shift necessary
end

v = kron(xvec,zvec);
w = kron(T,v);

u1 = 2*real(uspec1.'*w)+y;       % u1+y adds the laminar flow to u
u2 = 2*real(uspec2.'*w);
u3 = 2*real(uspec3.'*w);
u = [u1; u2; u3];   

v1a = kron(xvec2,zvec);
v1b = kron(T,v1a);
v2b = kron(derT,v);
v3a = kron(xvec,zvec2);
v3b = kron(T,v3a);

T11 = 2*real(uspec1.'*v1b);
T12 = 2*real(uspec1.'*v2b)+1 ;   % laminar addition
T13 = 2*real(uspec1.'*v3b);
T21 = 2*real(uspec2.'*v1b);
T22 = 2*real(uspec2.'*v2b);
T23 = 2*real(uspec2.'*v3b); 
T31 = 2*real(uspec3.'*v1b);
T32 = 2*real(uspec3.'*v2b);
T33 = 2*real(uspec3.'*v3b);
% velocity gradient tensor
T = [T11 T12 T13;T21 T22 T23; T31 T32 T33]
trace = T11 + T22 + T33
[evecs evals] = eig(T)



