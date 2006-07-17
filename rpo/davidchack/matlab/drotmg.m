function [g, jac] = drotmg(x0)

  global L M c P
  f = 8.0;  drm_init;
  
  g = x0;  
  if nargout > 1, jac = eye(4); end,
  for ip = 1:P,
    if nargout > 1, [g, jac1] = drotm(g);  jac = jac1*jac;
    else  g = drotm(g);  end,
  end,
  g = mod((g - x0+0.5*pi),2*pi)-0.5*pi;
  if nargout > 1, jac = jac - eye(4);  end,
return;
