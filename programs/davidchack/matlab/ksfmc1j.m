function dg = ksfmc1j(t,a)
%  function dg = ksfmc1j(t,a)

  persistent njeval
  
  N = 32;  nu = 0.015;  scale = 0.005;
  [f, df] = ksfmj(nu, a(1:N));  %nf = f./norm(f);
  b0 = [a(1:N) eye(N)];
  [b,s] = flowdp8('ksfmt', [N*(N+1); nu], b0(:), 0, [0 .5 1]*a(N+1), 2e-4);
  g = b(1:N,end) - a(1:N);
  ft = ksfm(nu, b(1:32,end));
  dg = reshape(b(N+1:end,end),N,N) - eye(N);
  dg = [dg ft; -scale.*[f'*dg+g'*df f'*ft]];
  
  if isempty(njeval), njeval = 0; end,  njeval = njeval + 1;  disp(njeval);
return;