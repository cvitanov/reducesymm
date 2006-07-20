function [tt, aa, da] = kscvfmj(nu, tend, a0, h, np)
% Solution of Kuramoto-Sivashinsky equation by ETDRK4 scheme
%
% KS equation in the form used by Cvitanovic:
%  w_s = (w^2)_y - w_yy - nu*w_yyyy, periodic BCs on [-pi, pi], odd solutions: w(-y,s) = -w(y,s)
%     In terms of Fourier modes: w(y,s) = i*sum(a_k(s)*exp(i*k*y)) with a_(-k) = -a_k, a_0 = 0;
%
% Transformation into form used by Kassam and Trefethen:
%
%  u_t = -u*u_x - u_xx - nu*u_xxxx, periodic BCs on [0, d], odd solutions: u(d-x,t) = -u(x,t)
%
%  y = (x/d - 1/2)*pi;  x = (y/pi + 1)*d/2
%
%  w_s = (w^2)_x*(d/2pi) - w_xx*(d/2pi)^2 - nu*w_xxxx*(d/2pi)^4
%
%  s = t*(2pi/d)^2;  dt/ds = (d/2pi)^2;   nu = (2pi/d)^2;  s = nu*t;
%
%  w_t = (w^2)_x/(d/2pi) - w_xx - w_xxxx
%
%  w_t = w*w_x*(4pi/d) - w_xx - w_xxxx
% 
%  (4pi/d)*w = -u;    v = fft(u);  v = -i*(4pi/d*N)*a;
% 
%   h - stepsize,  ns - number of steps,
%   np - output every np-th step (np = 0 - output only the final value)
%
%  Computation is based on v = fft(u), so linear term is diagonal.
%  Adapted from: AK Kassam and LN Trefethen, SISC 2005

%  clear all; load ksfm064upo;  a0 = a0(1:end-1);  h = 0.002;  np = 1;

  d = length(a0);
  N = 2*d+2;         % N should be preferably power of 2 (for faster FFT)
  v = -1i*[zeros(1,d+1); [a0 eye(d)];  zeros(1,d+1); -flipud([a0 eye(d)])]*N;
  k = [0:d 0 -d:-1]';                     % wave numbers
  L = k.^2 - nu.*k.^4;                    % Fourier multipliers
  
  E = exp(h*L);  E2 = exp(h*L/2);
  M = 16;                                 % no. of points for complex means
  r = exp(1i*pi*((1:M)-.5)/M);            % roots of unity
  LR = h*L(:,ones(M,1)) + r(ones(N,1),:);
  Q = h*real(mean((exp(LR/2)-1)./LR ,2));
  f1 = h*real(mean((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3 ,2));
  f2 = h*real(mean((2+LR+exp(LR).*(-2+LR))./LR.^3 ,2));
  f3 = h*real(mean((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3 ,2));
  
  E =  repmat(E,1,d+1);   E2 = repmat(E2,1,d+1);  Q = repmat(Q,1,d+1);
  f1 = repmat(f1,1,d+1);  f2 = repmat(f2,1,d+1);  f3 = repmat(f3,1,d+1);
  
  aa = a0;  tt = 0;  g = [-1i*k repmat(-2i*k,1,d)];
  ns = ceil(tend/h);%  ns = 2; 
  for n = 1:ns-1,
    t = n*h;                   rif = real(ifft(v));   Nv = g.*fft(repmat(rif(:,1),1,d+1).*rif);
    a = E2.*v + Q.*Nv;         rif = real(ifft(a));   Na = g.*fft(repmat(rif(:,1),1,d+1).*rif);
    b = E2.*v + Q.*Na;         rif = real(ifft(b));   Nb = g.*fft(repmat(rif(:,1),1,d+1).*rif);
    c = E2.*a + Q.*(2*Nb-Nv);  rif = real(ifft(c));   Nc = g.*fft(repmat(rif(:,1),1,d+1).*rif);
    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
    v = 0.5*(v(2:d+1,:)-flipud(v(d+3:end,:)));  v = [zeros(1,d+1); v; zeros(1,d+1); -flipud(v)];   % For odd solutions
    if np > 0 & mod(n,np)==0,  aa = [aa, -imag(v(2:d+1,1))/N];  tt = [tt, t];  end,
  end
  s = (tend-t)./h;  y0 = -imag(v(2:d+1,:))/N;
  
  rif = real(ifft(v));   Nv = g.*fft(repmat(rif(:,1),1,d+1).*rif);   % The last step
  a = E2.*v + Q.*Nv;         rif = real(ifft(a));   Na = g.*fft(repmat(rif(:,1),1,d+1).*rif);
  b = E2.*v + Q.*Na;         rif = real(ifft(b));   Nb = g.*fft(repmat(rif(:,1),1,d+1).*rif);
  c = E2.*a + Q.*(2*Nb-Nv);  rif = real(ifft(c));   Nc = g.*fft(repmat(rif(:,1),1,d+1).*rif);
  v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
  
  switch 2, 
  case 1,  %%% 1st order interpolation: y(s) = y0 + (y1-y0)*s for s in [0, 1],  dy/dt = y'(s)/h = (y1-y0)/h;
    y1 = -imag(v(2:d+1))/N;  ap = (y1-y0)/h;  ys = y0 + h*ap*s;
  case 2,  %%% 3rd order interpolation based on derivatives at end points: 
           %%% y(s) = y0 + h*y0p*s + (3*(y1-y0)-h*(y1p-2*y0p))*s^2 + (h*(y1p+y0p)-2*(y1-y0))*s^3;
           %%% dy/dt = y0p + 2*(3*(y1-y0)/h-(y1p-2*y0p))*s + 3*((y1p+y0p)-2*(y1-y0)/h)*s^2;
    y0p = repmat(L(2:d+1),1,d+1).*y0 - imag(Nv(2:d+1,:))/N;
    y1 = -imag(v(2:d+1,:))/N;  rif = real(ifft(v));   Nv = g.*fft(repmat(rif(:,1),1,d+1).*rif);
    y1p = repmat(L(2:d+1),1,d+1).*y1 - imag(Nv(2:d+1,:))/N;

    dy = (y1 - y0)./h;
    ys = y0 + h.*(y0p.*s + (3.*dy-y1p-2.*y0p).*s^2 + (y1p+y0p-2.*dy).*s^3);
  end,
  da = ys(:,2:d+1);
  
  if np > 0, aa = [aa, ys(:,1)];  tt = [tt, tend];  else aa = ys(:,1);  tt = tend; end,
