function g = ksfmf4(t, a, d, h, C, alpha);
%  function g = ksfmf4(t, a, d, h, C, alpha)  -  Associated flow for the full KSFM
%  4 - project away motion in the u_x direction

global TT A F FT NFEVAL
  
  [TT, A, F, FT] = ksfmedt(d, a(end), a(1:end-1), h, 2);

  k = 1:(length(a)-1)/2;
  b = reshape(a(1:end-1),2,[]);  b = reshape([b(2,:).*k; -b(1,:).*k],[],1);
  b = b./norm(b);  P = eye(length(a)-1) - b*b';
  
  g = P*(A(:,end) - a(1:end-1));
  g = [C*g; -alpha.*(FT'*g)];
  
  NFEVAL = NFEVAL + 1;
return;
