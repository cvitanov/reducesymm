function dg = ksfmc4j(t,a)
%  function dg = ksfmc4j(t,a)

  global C
  persistent njeval
  
  N = 32;  nu = 0.015;  alpha = 0.01;  h = 2e-4;
  [b, jac] = ksfmflowmapj(nu, a(N+1), a(1:N), h);
  g = C*(b - a(1:N));  dg = jac - eye(N);
  [ft, df] = ksfmj(nu, b);
  dg = [C*dg C*ft; -alpha.*[g'*df*jac+ft'*dg  g'*df*ft+ft'*ft]];
  
  if isempty(njeval), njeval = 0; end,  njeval = njeval + 1;%  disp(njeval);
return;