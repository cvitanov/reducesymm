function [f, df] = ksfd(t, u, L, sig);
%  Kuramoto-Sivashinsky equation
%                 u_t = -u*u_x - u_xx - u_xxxx,  x \in [0, L]
%    expressed in finite differences 
%    with symmetric (sig = 1) or antisymmetric (sig = -1) boundary conditions

  N = length(u);  dx = L./(N+1);
  y = [sig*u(1); 0; u(:); 0; sig*u(N)];
  c0 = 0.5/dx;  c3 = dx.^(-4);  c1 = 2.0./dx.^2-6.0*c3;  c2 = -1.0/dx.^2+4.0*c3;
  
  f = (c1 + c0.*(y(2:N+1)-y(4:N+3))).*y(3:N+2) + c2.*(y(2:N+1)+y(4:N+3)) - c3*(y(5:N+4)+y(1:N));

  if nargout > 1,
    df = zeros(N);
    df(1:(N+1):(N^2))=c1+c0.*(y(2:N+1)-y(4:N+3));  df(1,1)=df(1,1)-sig*c3;  df(N,N)=df(N,N)-sig*c3;
    df((N+1):(N+1):(N^2)) = -c0.*y(3:N+1) + c2;  df(2:(N+1):(N^2)) = c0.*y(4:N+2) + c2;
    df((2*N+1):(N+1):(N^2)) = -c3;  df(3:(N+1):(N^2-N)) = -c3;
  end,
  