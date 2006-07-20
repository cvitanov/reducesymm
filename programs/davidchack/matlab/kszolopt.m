function [g,dg] = kszolopt(x);
% function [g,dg] = kszolopt(x);
 
  global NFEVAL U0 UT TT
  N = 99;  dx = 0.5;  sig = 1;  h = 0.02;

  U0 = x(1:N);
  if nargout == 1,
    [UT, TT] = flowdp8('ks01', [N dx sig], U0, 0.0, x(N+1), h);
    g = UT - U0;  NFEVAL = NFEVAL + 1;
  else,
    dg = eye(N);  udu = [U0 dg];  udu = udu(:);  
    par = [N*(N+1) dx sig N];
    [udut,TT] = flowdp8('ks01t', par, udu, 0.0, x(N+1), h);
    UT = udut(1:N);  g = UT - U0;  ft = kszolflow(UT, dx, sig);
    dg = [reshape(udut(N+1:end),N,N)-eye(N) ft];
    NFEVAL = NFEVAL + N + 1;
  end