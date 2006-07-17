function [f, df] = ksfft(nu,a)
%KSFFT   Kuramoto-Sivashinski equation in Fourier modes for odd solutions: da/dt = f(a)
%   [f, df] = ksfft(nu,a) - using FFT for calculation
%   df = df/da - Jacobian

  d = length(a);  k = (1:d)';
  f = -(d+1).*real(fft(real(ifft(2i*[0; a; 0; -flipud(a)])).^2));
  f = k.^2.*(1-nu.*k.^2).*a(:) - k.*f(2:d+1);
  
  if nargout > 1,
    df = diag(k.^2.*(1-nu.*k.^2));
    for ik = k',  for ij = k',
      if ik > ij,  df(ik,ij) = df(ik,ij) - 2.*ik.*a(ik-ij);  end,
      if ik+ij <= d,  df(ik,ij) = df(ik,ij) + 2.*ik.*a(ik+ij);  end,
      if ij > ik,  df(ik,ij) = df(ik,ij) + 2.*ik.*a(ij-ik);  end, 
    end,  end,
  end,