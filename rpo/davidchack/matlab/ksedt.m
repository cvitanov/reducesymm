function [tt, uu] = ksedt(d, u, h, ns, np)
% Solution of Kuramoto-Sivashinsky equation by ETDRK4 scheme
%
% u_t = -u*u_x - u_xx - u_xxxx, periodic BCs on [0, d]
%   h - stepsize,  ns - number of steps,
%   np - output every np-th step (np = 0 - output only the final value)
%
%  Computation is based on v = fft(u), so linear term is diagonal.
%  Adapted from: AK Kassam and LN Trefethen, SISC 2005

  N = length(u);    % N should be even (preferably power of 2)
  v = fft(u);
  
  k = (2.*pi./d).*[0:N/2-1 0 -N/2+1:-1]'; % wave numbers
  L = k.^2 - k.^4;                        % Fourier multipliers
  E = exp(h*L);  E2 = exp(h*L/2);
  M = 16;                                 % no. of points for complex means
  r = exp(1i*pi*((1:M)-.5)/M);            % roots of unity
  LR = h*L(:,ones(M,1)) + r(ones(N,1),:);
  Q = h*real(mean((exp(LR/2)-1)./LR ,2));
  f1 = h*real(mean((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3 ,2));
  f2 = h*real(mean((2+LR+exp(LR).*(-2+LR))./LR.^3 ,2));
  f3 = h*real(mean((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3 ,2));

  uu = u;  tt = 0;  g = -0.5i*k;
  for n = 1:ns,
    t = n*h;                   Nv = g.*fft(real(ifft(v)).^2);
    a = E2.*v + Q.*Nv;         Na = g.*fft(real(ifft(a)).^2);
    b = E2.*v + Q.*Na;         Nb = g.*fft(real(ifft(b)).^2);
    c = E2.*a + Q.*(2*Nb-Nv);  Nc = g.*fft(real(ifft(c)).^2);
    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
    if np > 0 & mod(n,np)==0,
      u = real(ifft(v));  uu = [uu,u];  tt = [tt,t];
    elseif n == ns,
      uu = v;  tt = t;
    end
  end
