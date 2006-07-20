function g = kszolstb(s, x);
% function g = kszolstb(s, x);
  
  global C
  global NFEVAL U0 UT TT
  N = 50;  dx = 0.5;  sig = -1;  h = 0.02;  alpha = 1;

  U0 = x(1:N); 
  [UT, TT] = flowdp8('ks01', [N dx sig], U0, 0.0, x(N+1), h);
  ft = kszolflow(UT, dx, sig);

  g = UT - U0;  g = [C*g; -alpha.*(ft'*g)];

  NFEVAL = NFEVAL + 1;