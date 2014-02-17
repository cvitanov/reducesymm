%Adaptation of kursiv.m from Kassam ^ Trefethen (2005)

%Spatial grid:
d = 22
N = 32;
x = d*(1:N)'/N;
%Initial condition:
%u = cos(4*pi*x/d).*(1+sin(4*pi*x/d));
u = randn(32,1);
v = fft(u);

%ETDRK4 scalars:
h = 0.1;                        % time step
k = (2*pi./d)*[0:N/2-1 0 -N/2+1:-1]';  % wave numbers
%Linear term:
L = k.^2 - k.^4;                % Fourier multipliers

E = exp(h*L); E2 = exp(h*L/2);
M = 16;                         % no. of points for complex means
r = exp(1i*pi*((1:M)-.5)/M);    % roots of unity
LR = h*L(:,ones(M,1)) + r(ones(N,1), :);
Q  = h*real(mean(           (exp(LR/2) - 1)./LR              ,2));
f1 = h*real(mean(   (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3   ,2));
f2 = h*real(mean(       (2+LR+exp(LR).*(-2+LR))./LR.^3       ,2));
f3 = h*real(mean(   (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3   ,2));

%Time-stepping loop:
uu = u; tt = 0; vv = v;
tmax = 400; nmax = round(tmax/h); nplt = floor((tmax/1)/h);

%Nonlinear term:
g = -0.5i*k;
N = @(x) g.*fft(real(ifft(x)).^2);

for n = 1:nmax
    t = n*h;
    Nv = N(v);
    a = E2.*v + Q.*Nv;
    Na = N(a);
    b = E2.*v + Q.*Na;
    Nb = N(b);
    c = E2.*a + Q.*(2*Nb-Nv);
    Nc = N(c);
    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
    %if mod(n,nplt)==0
        u=real(ifft(v));
        uu = [uu, u]; tt = [tt, t]; vv = [vv, v];
    %end
end

figure(1)
plot3(real(vv(2,:)), imag(vv(2,:)), real(vv(3,:)))

#figure(2)
#surf(tt,x,uu)
#view([-90,90])
#shading('interp')
