function g = ksfmf6(t, a, d, h, C, alpha);
%  function g = ksfmf6(t, a, d, h, C, alpha)  -  Associated flow for the full KSFM
%  6 - determine integration time from dg/dT = 0 (local minimum)

global TT A F FT NFEVAL
  
  [TT, A, F, FT] = ksfmedt(d, a(end), a(1:end-1), h, 2);
    
  g = A(:,end) - a(1:end-1);
  g = [C*g; -alpha.*(FT'*g)];
  
  NFEVAL = NFEVAL + 1;
return;
