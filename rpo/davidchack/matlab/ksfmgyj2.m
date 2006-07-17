function [g, jac] = ksfmgyj2(y)
%  function [g, jac] = ksfmgyj2(y)
%     jac = dg/dy - optional

  N = 32;  nu = 0.015;  alpha = 0.005;  h = 2e-4;
  if nargout == 1,
    f = ksfm(nu, y(1:N));
    b = ksfmflowmapj(nu, y(N+1), y(1:N), h);
    g = b-y(1:N);  g = [g; -alpha.*(f'*g)];
  else
    [f, df] = ksfmj(nu, y(1:N));
    [b, jac] = ksfmflowmapj(nu, y(N+1), y(1:N), h);
    g = b-y(1:N);  ft = ksfm(nu, b);
    jac = [jac-eye(N) ft; -alpha.*[f'*jac+g'*df f'*ft]];
    g = [g; -alpha.*(f'*g)];
  end,
return;