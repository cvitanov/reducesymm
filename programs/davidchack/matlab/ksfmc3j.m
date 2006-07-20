function dg = ksfmc3j(t,a)
%  function dg = ksfmc3j(t,a)

  global C
  persistent njeval
  
  N = 32;  nu = 0.015;  alpha = 1.0;  h = 2e-4;
  [f, df] = ksfmj(nu, a(1:N));  %nf = f./norm(f);
  [b, dg] = ksfmflowmapj(nu, a(N+1), a(1:N), h);
  g = b - a(1:N);  dg = dg - eye(N);  ft = ksfmj(nu, b);
  dg = [C*dg C*ft; -alpha.*[f'*dg+g'*df f'*ft]];
  
  if isempty(njeval), njeval = 0; end,  njeval = njeval + 1;%  disp(njeval);
return;