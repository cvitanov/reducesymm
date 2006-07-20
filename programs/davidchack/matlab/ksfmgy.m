function g = ksfmgy(y)
%  function g = ksfmgy(y)

  N = 32;  nu = 0.015;  alpha = 0.005;
  f = ksfm(nu, y(1:N));
  [b,s] = flowdp8('ksfm', [N; nu], y(1:N), 0, [0 .5 1].*y(N+1), 2e-4);
  
  g = b(1:N,end)-y(1:N);  g = [g; -alpha.*(f'*g)];

return;