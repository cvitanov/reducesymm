function f = ksfm2(a, L)
%  function f = ksfm2(a, L) - a is a matrix
%
%KSFM   Kuramoto-Sivashinski equation in Fourier modes
%       u_t = -u*u_x - u_xx - u_xxxx, periodic BCs on [0, L]  
%  FM:  (u_k)_t = [(qk)^2-(qk)^4]*u_k + 0.5i*k*IFFT(FFT(u_k)^2)
%        q = 2pi/L;  
%  NOTE: Matlab defines FFT as inverse in most texts (e.g. Numerical Recipes)
  
  n = size(a,1); nt = size(a,2);
  ba = [zeros(1,nt); a(1:2:n-1,:)+1i*a(2:2:n,:); ...
        zeros(1,nt); a(n-1:-2:1,:)-1i*a(n:-2:2,:)];
  k = (2*pi/L).*[0:n/2, 0, -n/2:-1]';  ua = real(ifft(ba));
  bf = repmat(k.^2.*(1-k.^2),1,nt).*ba + 0.5i*(n+2)*repmat(k,1,nt).*fft(ua.^2);
  
  f = zeros(size(a));
  f(1:2:n-1,:) = real(bf(2:n/2+1,:));  f(2:2:n,:) = imag(bf(2:n/2+1,:));
  