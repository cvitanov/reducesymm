function [tt, aa] = ksfmmint(d, a0, tend, h, np)
% Solution of Kuramoto-Sivashinsky equation by ETDRK4 scheme
% ksfmmint - integrate till a local minimum of |aa(tt)-aa(0)| near tend
%
% u_t = -u*u_x - u_xx - u_xxxx, periodic BCs on [0, d]
%   h - stepsize
%   np - output every np-th step (np = 0 - output only the final value)
%
%  Computation is based on v = fft(u), so linear term is diagonal.
%  Adapted from: AK Kassam and LN Trefethen, SISC 2005

  if nargin < 5, np = 0; end
  N = length(a0)+2;  Nh = N/2;  % N should be even (preferably power of 2)
  v = [0; a0(1:2:end-1)+1i*a0(2:2:end); 0; a0(end-1:-2:1)-1i*a0(end:-2:2)];
  
  k = (2.*pi./d).*[0:Nh-1 0 -Nh+1:-1]';     % wave numbers
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
  g = 0.5i*k;  ns = ceil(tend/h);  
  tt = (0:ns)'.*h;  aa = zeros(N-2,ns+1);  aa(:,1) = a0;
  if nargout > 2,  Nv = g.*ifft(real(fft(v)).^2);
    vp = L(2:Nh).*v(2:Nh) + Nv(2:Nh);
    a0p = reshape([real(vp) imag(vp)]',[],1); end
  for n = 1:ns,
    t = n*h;                   Nv = g.*ifft(real(fft(v)).^2);
    a = E2.*v + Q.*Nv;         Na = g.*ifft(real(fft(a)).^2);
    b = E2.*v + Q.*Na;         Nb = g.*ifft(real(fft(b)).^2);
    c = E2.*a + Q.*(2*Nb-Nv);  Nc = g.*ifft(real(fft(c)).^2);
    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
    aa(1:2:end,n+1) = real(v(2:Nh));  aa(2:2:end,n+1) = imag(v(2:Nh));
  end
  
%% Check if tmin < tend+
  y1 = aa(:,end);  y1p = ksfm(0, y1, d);
  if y1p'*(y1-a0) >= 0,  % tmin < tend+: backtrack to find tmin
    for ii = size(aa,2):-1:2,
      y0 = aa(:,ii);  y0p = ksfm(0, y0, d);
      if y0p'*(y0-a0) < 0,  break;  else, y1 = y0; y1p = y0p; end, end
    aa = aa(:,1:ii);  tt = tt(1:ii);
  else, n = ns;          % tmin > tend+: take more steps to reach tmin
    while y1p'*(y1-a0) < 0,
      n = n + 1;  t = n*h;       Nv = g.*ifft(real(fft(v)).^2);
      a = E2.*v + Q.*Nv;         Na = g.*ifft(real(fft(a)).^2);
      b = E2.*v + Q.*Na;         Nb = g.*ifft(real(fft(b)).^2);
      c = E2.*a + Q.*(2*Nb-Nv);  Nc = g.*ifft(real(fft(c)).^2);
      v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
      y0 = y1;  y0p = y1p;  y1(1:2:end) = real(v(2:Nh));
      y1(2:2:end) = imag(v(2:Nh));  y1p = ksfm(0, y1, d);
      if n > ns, aa = [aa y0];  tt = [tt; t-h]; end, end, end
  
%% Find tmin using 3rd order interpolant (based on derivatives at end points) 
  dy = (y1 - y0)./h;  dya = (y0 - a0)./h;
  a2 = 3.*dy-y1p-2.*y0p;  a3 = y1p+y0p-2.*dy;
  cc = zeros(6,1);  cc(1) = y0p'*dya;
  cc(2) = 2.*a2'*dya + y0p'*y0p;   cc(3) = 3.*(a3'*dya + a2'*y0p);
  cc(4) = 2.*a2'*a2 + 4.*a3'*y0p;  cc(5) = 5.*a3'*a2;  cc(6) = 3.*a3'*a3;
  s0 = 1./(1-(y1p'*(y1-a0)./h)./cc(1)); % Initial guess from linear interpolant
  sz = fzero(@(s)polyval(flipud(cc),s),s0);
  ys = y0 + h.*(y0p.*sz + a2.*sz^2 + a3.*sz^3);
%  if nargout > 3, ap = y0p + 2.*(3.*dy-y1p-2.*y0p).*s + 3.*(y1p+y0p-2.*dy).*s^2; end 
  if np > 0, aa = [aa(:,1:np:end), ys];  tt = [tt(1:np:end); tt(end)+h.*sz]; 
  else aa = ys;  tt = tt(end)+h.*s; end
  