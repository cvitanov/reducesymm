function g = ksfmfm(t, a, d, h, C)
% function g = ksfmfm(t, a, d, h, C)
%   Associated flow for the full KSFM with minimum of |g| to determine tend

global TT A F FT NFEVAL
  
  ttry = TT(end); 
  [TT, A] = ksfmmint(d, a, ttry, h, 2);  g = C*(A(:,end) - a);
  
  NFEVAL = NFEVAL + 1;
