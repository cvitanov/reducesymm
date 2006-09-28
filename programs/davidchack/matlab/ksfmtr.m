function [f, df] = ksfmtr(a, L)
%  Traveling wave solutions in KSFM
%KSFM   Kuramoto-Sivashinski equation in Fourier modes
%       u_t = -u*u_x - u_xx - u_xxxx, periodic BCs on [0, L]  
%  FM:  (u_k)_t = [(qk)^2-(qk)^4]*u_k + 0.5i*k*IFFT(FFT(u_k)^2)
%        q = 2pi/L;  
%  NOTE: Matlab defines FFT as inverse in most texts (e.g. Numerical Recipes)
  
  n = length(a)-1;  c = a(end);  k = (2*pi/L)*(1:n/2)';
  if nargout > 1,  [f, df] = ksfm(0, a(1:n), L);
  else f = ksfm(0, a(1:n), L); end,
  f(1:2:n) = f(1:2:n) + c*k.*a(2:2:n); f(2:2:n) = f(2:2:n) - c*k.*a(1:2:n);
  f = [f; f(1:2:n)'*(k.*a(2:2:n)) - f(2:2:n)'*(k.*a(1:2:n))];

  if nargout > 1,
%    disp(n); disp(size(df)); disp(size(df(2:2*n+2:end))); disp(size(c*k));
    fc = (k.*a(2:2:n))'*df(1:2:n,1:n)-(k.*a(1:2:n))'*df(2:2:n,1:n);
    df(2:2*n+2:end) = df(2:2*n+2:end) - c*k';
    df(n+1:2*n+2:end) = df(n+1:2*n+2:end) + c*k';
    df = [df zeros(n,1); zeros(1,n+1)];
    df(1:2:n,n+1) = k.*a(2:2:n);  df(2:2:n,n+1) = -k.*a(1:2:n);
%    df(n+1,1:n) = (k.*a(2:2:n))'*df(1:2:n,1:n)-(k.*a(1:2:n))'*df(2:2:n,1:n);
%    disp(df(n+1,1:n));  disp((k.*f(2:2:n))');
    df(n+1,1:2:n) = fc(1:2:n) - (k.*f(2:2:n))' + c*(k.^2.*a(1:2:n))';
    df(n+1,2:2:n) = fc(2:2:n) + (k.*f(1:2:n))' + c*(k.^2.*a(2:2:n))';
%    fc = f(1:n);  fc(2:2:n) = -fc(2:2:n);
%    df(n+1,1:n) = fc'*df(1:n,1:n);
    df(n+1,n+1) = (k.^2)'*(a(1:2:n).^2 + a(2:2:n).^2);
  end