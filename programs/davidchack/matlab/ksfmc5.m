function g = ksfmc5(t, a);
%  function g = ksfmc5(t, a);

global C 
global TT A F FT NFEVAL
  
%  N = 31;  nu = (2*pi/38.5)^2;  h = 0.002; alpha = 0.01; 
  N = 63;  nu = (2*pi/54.0)^2;  h = 0.002; alpha = 0.01; 
  [TT, A, F, FT] = kscvfm(nu, a(N+1), a(1:N), h, 3);
    
  g = A(:,end) - a(1:N);
  g = [C*g; -alpha.*(FT'*g)];
  
  NFEVAL = NFEVAL + 1;
return;
