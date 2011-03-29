function xz = ksfm2zeros(a, L, flg)

  if nargin < 3, flg = 0; end
  switch flg,
  case 0, [x, u] = ksfm2real(a, L);
  case 1, [x, u, ux] = ksfm2real(a, L); u = ux;
  case 2, [x, u, ux, uxx] = ksfm2real(a, L); u = uxx; end
  ix = find((u(1:end-1)<0) == (u(2:end)>=0));
  
  xz = zeros(size(ix));
  for ii = 1:length(ix),
    xz(ii) = fzero(@(x)ksfm2anyx(a,L,x,flg), x(ix(ii)+[0 1]));
%    [xz(ii), fval] = fzero(@(x)ksfm2anyx(a,L,x,flg), x(ix(ii)+[0 1]));
%    disp([xz(ii) fval]); pause;
  end