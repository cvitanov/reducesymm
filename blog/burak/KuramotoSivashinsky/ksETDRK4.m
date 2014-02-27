%Adaptation of kursiv.m from Kassam ^ Trefethen (2005)

%Spatial grid:
N = 128;
x = 32*pi*(1:N)'/N;
%Initial condition:
u = cos(x/16).*(1+sin(x/16));
v = fft(u);
v = randn(128,1) + 1i*randn(128,1);
v(1) = 0;
v(2) = 0;
v(64) = 0;
v(end) = 0;

%ETDRK4 scalars:
h = 1/4;                        % time step
k = [0:N/2-1 0 -N/2+1:-1]'/16;  % wave numbers
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
tmax = 150; nmax = round(tmax/h); nplt = 1; %floor((tmax/100)/h);
g = -0.5i*k;

%Nonlinear term:
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
    if mod(n,nplt)==0
        u=real(ifft(v));
        uu = [uu, u]; tt = [tt, t];; vv = [vv, v];
    end
end

figure(1)
plot(real(vv(2,:)))

figure(2)
surf(tt,x,uu)
view([-90,90])
shading('interp')
