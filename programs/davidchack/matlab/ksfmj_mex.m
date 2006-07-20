function da = ksfmj_mex(t,a)

  k = 1:length(a);  nu = 0.015;
  da = diag(k.^2.*(1-nu.*k.^2));