function g = ksfmf2(t, a, d, h, C, alpha);
%  function g = ksfmf2(t, a, d, h, C, alpha) - Associated flow for the full KSFM
%  2 - constrain A(1) = 0, so a = [A(2:end); tend]

global TT A F FT NFEVAL
  
  [TT, A, F, FT] = ksfmedt(d, a(end), [0; a(1:end-1)], h, 2);
    
  g = A(2:end,end) - a(1:end-1);
  g = [C*g; -alpha.*(FT(2:end)'*g)];
  
  NFEVAL = NFEVAL + 1;
return;
