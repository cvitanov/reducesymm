function f = ksfmtr(a, d)
%  Traveling wave solutions in KSFM
%KSFM   Kuramoto-Sivashinski equation in Fourier modes
%       u_t = -u*u_x - u_xx - u_xxxx, periodic BCs on [0, d]  
%  FM:  (u_k)_t = [(qk)^2-(qk)^4]*u_k + 0.5i*k*IFFT(FFT(u_k)^2)
%        q = 2pi/d;  
%  NOTE: Matlab defines FFT as inverse in most texts (e.g. Numerical Recipes)
  
  n = length(a) - 1;  c = a(end);
  ba = [0; a(1:2:n-1)+1i*a(2:2:n); 0; a(n-1:-2:1)-1i*a(n:-2:2)];
  k = (2*pi/d).*[0:n/2, 0, -n/2:-1]';  ua = real(fft(ba));
  bf = k.^2.*(1-k.^2).*ba + 0.5i*k.*ifft(ua.^2) - 1i*c*k.*ba;
  
  f = zeros(size(a));  
  f(1:2:n-1) = real(bf(2:n/2+1));
  f(2:2:n) = imag(bf(2:n/2+1));
  f(end) = real(conj(bf)'*(1i*k.*ba));
