% u_sum.m

function u = u_sum(tt,rr)

global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N uspec1 uspec2 uspec3

x = rr(1);
y = rr(2);
z = rr(3);

zvec = [1/2];
for mz = 1:Mz-1
    zvec = [zvec; exp(2*pi*i*mz*z/Lz)];
end

xvec = [];
for mx = 0:Mx-1
    if 0 <= mx & mx <= Mx/2
        kx = mx;
    else
        kx = mx - Mx;
    end      
    xvec = [xvec; exp(2*pi*i*kx*x/Lx)];
end

T = zeros(My,1);                    % should generalize to [a,b].
T(1) = 1;        % 0th cheby poly
T(2) = y;
for my = 3:My
    T(my) = 2*y*T(my-1) - T(my-2);
end

v = kron(xvec,zvec);
w = kron(T,v);

u1 = 2*real(uspec1.'*w);
u2 = 2*real(uspec2.'*w);
u3 = 2*real(uspec3.'*w);
u = [u1+y; u2; u3];         % u1+y adds the laminar flow to u
