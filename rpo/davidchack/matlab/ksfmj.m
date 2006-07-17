function [f, df] = ksfmj(nu,a)
%KSFMJ   Kuramoto-Sivashinski equation in Fourier modes for odd solutions: da/dt = f(a)
%   [f, df] = ksfmj(nu,a);
%   df = df/da - Jacobian

  d = length(a);  k = (1:d)';
  f = k.^2.*(1-nu.*k.^2).*a(:);
  for ik = k',
    f(ik) = f(ik) - ik.*(sum(a(1:ik-1).*a(ik-1:-1:1)) - 2.*sum(a(1:d-ik).*a(ik+1:d)));
  end,
  
%  ff = -(2*d+2).*real(fft(real(ifft(1i*[0; a; 0; -flipud(a)])).^2));
%  f = k.^2.*(1-nu.*k.^2).*a(:) - k.*ff(2:d+1);
  
  if nargout > 1, 
    df = diag(k.^2.*(1-nu.*k.^2));
    for ik = k',  for ij = k',
      if ik > ij,  df(ik,ij) = df(ik,ij) - 2.*ik.*a(ik-ij);  end,
      if ik+ij <= d,  df(ik,ij) = df(ik,ij) + 2.*ik.*a(ik+ij);  end,
      if ij > ik,  df(ik,ij) = df(ik,ij) + 2.*ik.*a(ij-ik);  end, 
    end,  end,
  end,