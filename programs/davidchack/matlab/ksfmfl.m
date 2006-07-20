function g = ksfmfl(t, a, d, h, C)
% function g = ksfmfl(t, a, d, h, C)
%   Associated flow for the full KSFM with minimum of |g| to determine tend
%   Add motion along the flow

global TT A F FT NFEVAL
  
  ttry = TT(end); 
  [TT, A] = ksfmmint(d, a, ttry, h, 2);  F = ksfm(0, a, d);  
  g = A(:,end) - a;  ee = 1.0*g'*g;  g = C*g + ee.*F;
  
  NFEVAL = NFEVAL + 1;
