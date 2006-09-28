function [f, df] = ksfm(t, a, L, C)
%KSFM   Kuramoto-Sivashinski equation in Fourier modes
%       u_t = -u*u_x - u_xx - u_xxxx, periodic BCs on [0, L]  
%  FM:  (u_k)_t = [(qk)^2-(qk)^4]*u_k + 0.5i*k*IFFT(FFT(u_k)^2)
%        q = 2pi/L;  
%  NOTE: Matlab defines FFT as inverse in most texts (e.g. Numerical Recipes)
  
  n = length(a);
  ba = [0; a(1:2:n-1)+1i*a(2:2:n); 0; a(n-1:-2:1)-1i*a(n:-2:2)];
  k = (2*pi/L).*[0:n/2, 0, -n/2:-1]';  ua = real(fft(ba));
  bf = k.^2.*(1-k.^2).*ba + 0.5i*k.*ifft(ua.^2);
  
  f = zeros(size(a));
  f(1:2:n-1) = real(bf(2:n/2+1));  f(2:2:n) = imag(bf(2:n/2+1));
  
  if nargin > 3,  f = C*f;
%    disp([t norm(f)]);  pause;
  end
  
  if nargout > 1,
    bf = zeros(n+2,n);
    bf(2:n/2+1,1:2:end-1) = eye(n/2);    bf(n/2+3:end,end-1:-2:1) = eye(n/2);
    bf(2:n/2+1,2:2:end) = 1i*eye(n/2);   bf(n/2+3:end,end:-2:2) = -1i*eye(n/2);
    
    bf = repmat(k.^2.*(1-k.^2),1,n).*bf + 1i*repmat(k,1,n).*ifft(repmat(ua,1,n).*fft(bf));
    
    df = zeros(n);
    df(1:2:n-1,:) = real(bf(2:n/2+1,:));  df(2:2:n,:) = imag(bf(2:n/2+1,:));

    if nargin > 3, df = C*df; end
  end
