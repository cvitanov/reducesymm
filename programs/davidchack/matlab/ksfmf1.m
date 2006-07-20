function g = ksfmf1(t, a, d, h, C, alpha);
%  function g = ksfmf1(t, a, d, h, C, alpha)  -  Associated flow for the full KSFM

global TT A F FT NFEVAL
  
  [TT, A, F, FT] = ksfmedt(d, a(end), a(1:end-1), h, 2);
    
  g = A(:,end) - a(1:end-1);
  g = [C*g; -alpha.*(FT'*g)];
  
  NFEVAL = NFEVAL + 1;
return;
