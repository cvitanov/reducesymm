function [u, du] = lj(r);

  u = 4.*(r.^(-12) - r.^(-6));
  if nargout > 1,  du = -48.*(r.^(-13) - 0.5.*r.^(-7)); end,

return;