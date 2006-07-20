function dadt = ksfm_mex(t,a)
%KSFM   Kuramoto-Sivashinski equation in Fourier modes for odd solutions

  nu = 0.015;  d = length(a);  k = (1:d)';
  dadt = k.^2.*(1-nu.*k.^2).*a(:);
  for ik = k',
    dadt(ik) = dadt(ik) - ik.*(sum(a(1:ik-1).*a(ik-1:-1:1)) - 2.*sum(a(1:d-ik).*a(ik+1:d)));
  end,
