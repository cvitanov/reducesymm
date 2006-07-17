function [g, jac] = ksfmgyj(y)
%  function [g, jac] = ksfmgyj(y)
%     jac = dg/dy - optional

  N = 32;  nu = 0.015;  alpha = 0.005;
  if nargout == 1,
    f = ksfm(nu, y(1:N));
    [b,s] = flowdp8('ksfm', [N; nu], y(1:N), 0, [0 .5 1].*y(N+1), 2e-4);
    g = b(1:N,end)-y(1:N);  g = [g; -alpha.*(f'*g)];
  else
    [f, df] = ksfmj(nu, y(1:N));
    b0 = [y(1:N) eye(N)];
    [b,s] = flowdp8('ksfmt', [N*(N+1); nu], b0(:), 0, [0 .5 1]*y(N+1), 2e-4);
    g = b(1:N,end) - y(1:N);
    ft = ksfm(nu, b(1:32,end));
    jac = reshape(b(N+1:end,end),N,N) - eye(N);
    jac = [jac ft; -alpha.*[f'*jac+g'*df f'*ft]];
    g = [g; -alpha.*(f'*g)];
  end,
return;