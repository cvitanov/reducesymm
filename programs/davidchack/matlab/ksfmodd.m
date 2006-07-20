function [f, df] = ksfmodd(nu,a)
%KSFFT   Kuramoto-Sivashinski equation in Fourier modes for odd solutions: da/dt = f(a)
%        a_k' = (k^2 - nu*k^4)*a_k - k*sum_m a_m*a_{k-m}
% 
%   [f, df] = kfmodd(nu,a) - using FFT for calculation
%   df = df/da - Jacobian

  d = length(a);  k = (1:d)';
  
  ff = real(ifft(1i*[0; a(:); 0; -flipud(a(:))]));
  
  f = -4.*(d+1).*real(fft(ff.^2));
  f = k.^2.*(1-nu.*k.^2).*a(:) - k.*f(2:d+1);

  
  if nargout > 1,
    
    df = real(ifft(1i*[zeros(1,d); eye(d); zeros(1,d); -flipud(eye(d))]));
    
    df = 8.*(d+1).*fft(diag(ff)*df);
    
    df = diag(k.^2.*(1-nu.*k.^2)) + repmat(k,1,d).*df(2:d+1,:);    
    
  end,