function [tt, aa, da] = ksfmjedt(d, tend, a0, h, np)
% Solution (with Jacobian) of Kuramoto-Sivashinsky equation by ETDRK4 scheme
%
% u_t = -u*u_x - u_xx - u_xxxx, periodic BCs on [0, d]
%   h - stepsize,  ns - number of steps,
%   np - output every np-th step (np = 0 - output only the final value)
%
%  Computation is based on v = fft(u), so linear term is diagonal.
%  Adapted from: AK Kassam and LN Trefethen, SISC 2005

  if nargin < 5, np = 0; end
  N = length(a0)+2;  Nh = N/2;  % N should be even (preferably power of 2)
  v = [zeros(1,N-1); [a0(1:2:end-1)+1i*a0(2:2:end) zeros(Nh-1,N-2)]; ...
       zeros(1,N-1); [a0(end-1:-2:1)-1i*a0(end:-2:2) zeros(Nh-1,N-2)]];
  v(N+2:2*N+1:end) = 1;  v(2*N+2:2*N+1:end) = 1i;
  v(2*N:2*N-1:end) = 1;  v(3*N:2*N-1:end) = -1i;
  
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

  E =  repmat(E,1,N-1);   E2 = repmat(E2,1,N-1);  Q = repmat(Q,1,N-1);
  f1 = repmat(f1,1,N-1);  f2 = repmat(f2,1,N-1);  f3 = repmat(f3,1,N-1);

  aa = a0;  tt = 0;  g = 0.5i*k;  ns = ceil(tend/h);
%  Nv = g.*ifft(real(fft(v(:,1))).^2);  vp = L(2:Nh).*v(2:Nh,1) + Nv(2:Nh);  
%  a0p = reshape([real(vp) imag(vp)]',[],1);
  g = [g repmat(2.*g,1,N-2)];
  for n = 1:ns,
    t = n*h;                   rfv = real(fft(v));   Nv = g.*ifft(repmat(rfv(:,1),1,N-1).*rfv);
    a = E2.*v + Q.*Nv;         rfv = real(fft(a));   Na = g.*ifft(repmat(rfv(:,1),1,N-1).*rfv);
    b = E2.*v + Q.*Na;         rfv = real(fft(b));   Nb = g.*ifft(repmat(rfv(:,1),1,N-1).*rfv);
    c = E2.*a + Q.*(2*Nb-Nv);  rfv = real(fft(c));   Nc = g.*ifft(repmat(rfv(:,1),1,N-1).*rfv);
    if n == ns, v0 = v; end,     % save penultimate 
    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
    if np > 0 & mod(n,np)==0 & n < ns,
      y1 = [real(v(2:Nh,1)) imag(v(2:Nh,1))]'; aa = [aa, y1(:)];  tt = [tt,t]; end
  end

%%% 3rd order interpolation based on derivatives at end points: 
%%% y(s) = y0 + h*y0p*s + (3*(y1-y0)-h*(y1p-2*y0p))*s^2 + (h*(y1p+y0p)-2*(y1-y0))*s^3;
%%% dy/dt = y0p + 2*(3*(y1-y0)/h-(y1p-2*y0p))*s + 3*((y1p+y0p)-2*(y1-y0)/h)*s^2;
  s = (tend-t)./h + 1;
  
  v0p = repmat(L(2:Nh),1,N-1).*v0(2:Nh,:) + Nv(2:Nh,:);  
  rfv = real(fft(v));   Nv = g.*ifft(repmat(rfv(:,1),1,N-1).*rfv);
  vp = repmat(L(2:Nh),1,N-1).*v(2:Nh,:) + Nv(2:Nh,:);
  y0 = zeros(N-2,N-1);  y0(1:2:end,:) = real(v0(2:Nh,:));  y0(2:2:end,:) = imag(v0(2:Nh,:));
  y0p = zeros(N-2,N-1);  y0p(1:2:end,:) = real(v0p);  y0p(2:2:end,:) = imag(v0p);
  y1 = zeros(N-2,N-1);  y1(1:2:end,:) = real(v(2:Nh,:));  y1(2:2:end,:) = imag(v(2:Nh,:));
  y1p = zeros(N-2,N-1);  y1p(1:2:end,:) = real(vp);  y1p(2:2:end,:) = imag(vp);
  
  dy = (y1 - y0)./h;
  ys = y0 + h.*(y0p.*s + (3.*dy-y1p-2.*y0p).*s^2 + (y1p+y0p-2.*dy).*s^3);
  da = ys(:,2:end);
  
%  ap = y0p(:,1) + 2.*(3.*dy(:,1)-y1p(:,1)-2.*y0p(:,1)).*s + 3.*(y1p(:,1)+y0p(:,1)-2.*dy(:,1)).*s^2;
  
  if np > 0, aa = [aa, ys(:,1)];  tt = [tt, tend]; else aa = ys(:,1);  tt = tend; end
  