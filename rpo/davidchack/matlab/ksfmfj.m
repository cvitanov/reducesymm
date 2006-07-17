function [f, df] = ksfmfj(t, a, nu, C)
%KSFMFJ   Kuramoto-Sivashinski equation in Fourier modes: da/dt = f(a)
%         f(a) = (k^2 - nu*k^4)*b(a) - 1i*k*FFT(IFFT(b(a)).^2)
%         where b(a) = [0; a(1:2:d-1)+1i*a(2:2:d); 0; a(d-1:-2:1)-1i*a(d:-2:2)],
%         k = [0:d/2; 0; -d/2:-1]';  d = length(a) (d = 30, 62, ... 2^n-2)
%   [f, df] = ksfmfj(t, a, nu, C);
%   df = df/da - Jacobian

%  nu = 0.02;  a = [1; .8; .5; -.3; .2; .1];

  d = length(a);
  ba = [0; a(1:2:d-1)+1i*a(2:2:d); 0; a(d-1:-2:1)-1i*a(d:-2:2)];
  k = [0:d/2, 0, -d/2:-1]';  ua = real(ifft(ba));
  bf = k.^2.*(1-nu.*k.^2).*ba - 1i*k.*fft(ua.^2).*(d+2);
  
  f = zeros(size(a));
  f(1:2:d-1) = real(bf(2:d/2+1));  f(2:2:d) = imag(bf(2:d/2+1));
  
  if nargin > 3,
    f = C*f;
    disp([t norm(f)]);  pause;
  end,
  
  if nargout > 1,
    bf = zeros(d+2,d);
    bf(2:d/2+1,1:2:end-1) = eye(d/2);    bf(d/2+3:end,end-1:-2:1) = eye(d/2);
    bf(2:d/2+1,2:2:end) = 1i*eye(d/2);   bf(d/2+3:end,end:-2:2) = -1i*eye(d/2);
    
    bf = repmat(k.^2.*(1-nu.*k.^2),1,d).*bf - 2i*repmat(k,1,d).*fft(repmat(ua,1,d).*ifft(bf)).*(d+2);
    
    df = zeros(d);
    df(1:2:d-1,:) = real(bf(2:d/2+1,:));  df(2:2:d,:) = imag(bf(2:d/2+1,:));
  end,
