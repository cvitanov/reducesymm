function [g, jac] = ksfmf5(t, a, d, h, C, alpha);
%  function g = ksfmf5(t, a, d, h, C, alpha)  -  Associated flow for the full KSFM
%  5 - for quasiperiodic orbits

global TT A F FT NFEVAL PHASE DD
  
  if nargout == 1,
    [TT, A, F, FT] = ksfmedt(d, a(end-1), a(1:end-2), h, 2);
  else
    [TT, A, F, FT, DA] = ksfmjedt(d, a(end-1), a(1:end-2), h, 2); end
  
  k = (1:length(a)/2-1)';  ek = exp(-2i*pi/d.*k.*a(end));
  v = (A(1:2:end-1,end)+1i*A(2:2:end,end)).*ek;
  
  g = a(1:end-2);  g(1:2:end) = real(v) - g(1:2:end);  
  g(2:2:end) = imag(v) - g(2:2:end); 
  
  vf = (FT(1:2:end-1)+1i*FT(2:2:end)).*ek;  ftp = zeros(size(g));  
  ftp(1:2:end-1) = real(vf);  ftp(2:2:end) = imag(vf);
  
  v = 1i*k.*v;  app = zeros(size(g));
  app(1:2:end-1) = real(v);  app(2:2:end) = imag(v);
  
  g = [C*g; -alpha(1).*(ftp'*g); alpha(2).*app'*g];
  
  NFEVAL = NFEVAL + 1;  PHASE = a(end);  DD = d;
return;
