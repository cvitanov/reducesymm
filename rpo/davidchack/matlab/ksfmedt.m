function [tt, aa, a0p, ap] = ksfmedt(d, tend, a0, h, np)
% Solution of Kuramoto-Sivashinsky equation by ETDRK4 scheme
%
% u_t = -u*u_x - u_xx - u_xxxx, periodic BCs on [0, d]
%   h - stepsize,  ns - number of steps,
%   np - output every np-th step (np = 0 - output only the final value)
%
%  Computation is based on v = fft(u), so linear term is diagonal.
%  Adapted from: AK Kassam and LN Trefethen, SISC 2005

  if nargin < 5, np = 0; end
  N = length(a0)+2;  Nh = N/2;  % N should be even (preferably power of 2)
  v = [0; a0(1:2:end-1)+1i*a0(2:2:end); 0; a0(end-1:-2:1)-1i*a0(end:-2:2)];
  
  k = (2.*pi./d).*[0:Nh-1 0 -Nh+1:-1]';   % wave numbers
  L = k.^2 - k.^4;                          % Fourier multipliers
  E = exp(h*L);  E2 = exp(h*L/2);
  M = 16;                                   % no. of points for complex means
  r = exp(1i*pi*((1:M)-.5)/M);              % roots of unity
  LR = h*L(:,ones(M,1)) + r(ones(N,1),:);
  Q = h*real(mean((exp(LR/2)-1)./LR ,2));
  f1 = h*real(mean((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3 ,2));
  f2 = h*real(mean((2+LR+exp(LR).*(-2+LR))./LR.^3 ,2));
  f3 = h*real(mean((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3 ,2));
% disp(tend);
  aa = a0;  tt = 0;  g = 0.5i*k;  ns = ceil(tend/h);
  if nargout > 2,  Nv = g.*ifft(real(fft(v)).^2);
    vp = L(2:Nh).*v(2:Nh) + Nv(2:Nh);  a0p = reshape([real(vp) imag(vp)]',[],1); end
  for n = 1:ns,
    t = n*h;                   Nv = g.*ifft(real(fft(v)).^2);
    a = E2.*v + Q.*Nv;         Na = g.*ifft(real(fft(a)).^2);
    b = E2.*v + Q.*Na;         Nb = g.*ifft(real(fft(b)).^2);
    c = E2.*a + Q.*(2*Nb-Nv);  Nc = g.*ifft(real(fft(c)).^2);
    if n == ns, v0 = v; end,     % save penultimate 
    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
    if np > 0 & mod(n,np)==0 & n < ns,
      y1 = [real(v(2:Nh)) imag(v(2:Nh))]'; aa = [aa, y1(:)];  tt = [tt,t]; end
  end

%%% 3rd order interpolation based on derivatives at end points: 
%%% y(s) = y0 + h*y0p*s + (3*(y1-y0)-h*(y1p-2*y0p))*s^2 + (h*(y1p+y0p)-2*(y1-y0))*s^3;
%%% dy/dt = y0p + 2*(3*(y1-y0)/h-(y1p-2*y0p))*s + 3*((y1p+y0p)-2*(y1-y0)/h)*s^2;
  s = (tend-t)./h + 1;
  
  v0p = L(2:Nh).*v0(2:Nh) + Nv(2:Nh);  Nv = g.*ifft(real(fft(v)).^2);
  vp = L(2:Nh).*v(2:Nh) + Nv(2:Nh);
  y0 = reshape([real(v0(2:Nh)) imag(v0(2:Nh))]',[],1);
  y0p = reshape([real(v0p) imag(v0p)]',[],1);
  y1 = reshape([real(v(2:Nh)) imag(v(2:Nh))]',[],1);
  y1p = reshape([real(vp) imag(vp)]',[],1);

  dy = (y1 - y0)./h;
  ys = y0 + h.*(y0p.*s + (3.*dy-y1p-2.*y0p).*s^2 + (y1p+y0p-2.*dy).*s^3);
  if nargout > 3, ap = y0p + 2.*(3.*dy-y1p-2.*y0p).*s + 3.*(y1p+y0p-2.*dy).*s^2; end
  
  if np > 0, aa = [aa, ys];  tt = [tt, tend]; else aa = ys;  tt = tend; end
  