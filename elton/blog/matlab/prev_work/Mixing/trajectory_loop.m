
Mx = 48; My = 35; Mz = 25; Nx = 48; Ny = 35; Nz = 48; Nd = 3; a = -1; b = 1;
Lx = 2*pi/1.14; Lz =  2*pi/2.5;

N = size(uNB_test,1);

uspec1 = zeros(N/3,1);
uspec1 = uNB_test(1:N/3,1) + i*uNB_test(1:N/3,2);

uspec2 = zeros(N/3,1);
uspec2 = uNB_test(N/3+1:2*N/3,1) + i*uNB_test(N/3+1:2*N/3,2);

uspec3 = zeros(N/3,1);
uspec3 = uNB_test(2*N/3+1:N,1) + i*uNB_test(2*N/3+1:N,2);

x(1) = pi/(1.14*8); y(1) = cos(5*pi/34); z(1) = pi/10;

dt = 0.1;
tf = 100;

tic
for jj = 1:1

zvec = [];
for mz = 0:Mz-1
    zvec = [zvec; exp(2*pi*i*mz*z(jj)/Lz)];
end

xvec = [];
for mx = 0:Mx-1
    if 0 <= mx & mx <= Mx/2
        kx = mx;
    else
        kx = mx - Mx;
    end      
    xvec = [xvec; exp(2*pi*i*kx*x(jj)/Lx)];
end

T = zeros(My,1);                    % should generalize to [a,b]
T(1) = 1;
T(2) = y(jj);
for my = 3:My
    T(my) = 2*y(jj)*T(my-1) - T(my-2);
end


offset = 1; u = 0 + 0*i;
for my = 1:My
    sum = 0 + 0*i;
    for mx = 1:Mx
        sum = sum + (zvec.'*uspec1(offset:Mz-1+offset))*xvec(mx);   
        offset = offset + Mz;
    end
    u = u + sum*T(my);
end
u1 = 2*real(u);

offset = 1; u = 0 + 0*i;
for my = 1:My
    sum = 0 + 0*i;
    for mx = 1:Mx
        sum = sum + (zvec.'*uspec2(offset:Mz-1+offset))*xvec(mx);   
        offset = offset + Mz;
    end
    u = u + sum*T(my);
end
u2 = 2*real(u);

offset = 1; u = 0 + 0*i;
for my = 1:My
    sum = 0 + 0*i;
    for mx = 1:Mx
        sum = sum + (zvec.'*uspec3(offset:Mz-1+offset))*xvec(mx);   
        offset = offset + Mz;
    end
    u = u + sum*T(my);
end
u3 = 2*real(u);

uu = [u1; u2; u3];


end 

toc









