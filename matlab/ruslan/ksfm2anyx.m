function [u, ux, uxx] = ksfm2anyx(a, L, x, flg)
%  [U, UX, UXX] = KSFM2ANYX(A, L, X)
%  Convert FM representation A to real space U(X) for any values of X
%  (i.e. not using FFT).
%  Also calculate U'(X) and U''(X) (optional).

  if nargin < 4, flg = 0; end
  n = length(a)/2;  nx = length(x);
  qk = (2*pi/L)*(1:n)';  qkx = qk*x(:)';
  coqkx = cos(qkx);  siqkx = sin(qkx);
  u = 2.*sum(repmat(a(1:2:end),1,nx).*coqkx + repmat(a(2:2:end),1,nx).*siqkx)';
  if nargout > 1 || flg == 1,
    ux = 2.*sum(-repmat(qk.*a(1:2:end),1,nx).*siqkx + repmat(qk.*a(2:2:end),1,nx).*coqkx)';  end
  if nargout > 2 || flg == 2,
    uxx = -2.*sum(repmat(qk.^2.*a(1:2:end),1,nx).*coqkx + repmat(qk.^2.*a(2:2:end),1,nx).*siqkx)';  end
  if flg == 1, u = ux; end
  if flg == 2, u = uxx; end
  