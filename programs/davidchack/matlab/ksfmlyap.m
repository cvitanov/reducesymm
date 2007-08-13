function [at, vt] = ksfmlyap(a0, v0, d, h, nstp)
% Solution of Kuramoto-Sivashinsky equation together with
% variational equation for Lyapunov exponents calculation
% by ETDRK4 scheme for a given number of steps
%
% u_t = -u*u_x - u_xx - u_xxxx, periodic BCs on [-d/2, d/2]
% 
% Input:
%   a0 - initial condition, v0 - initial condition for variational equation
%   d - length parameter
%   h - stepsize,  nstp - number of steps
%   np - output every np-th step (np = 0 - output only the final value)
% Output:
%   at - solution
%   vt - solution of the variational problem

%  Computation is based on v = fft(u), so linear term is diagonal.
%  Adapted from: AK Kassam and LN Trefethen, SISC 2005

  N = length(a0)+2;  Nh = N/2;  % N should be even (preferably power of 2)
  NN = size(v0,2);
  v = [zeros(1,NN+1); [a0(1:2:end-1)+1i*a0(2:2:end) v0(1:2:end-1,:)+1i*v0(2:2:end,:)];
       zeros(1,NN+1); [a0(end-1:-2:1)-1i*a0(end:-2:2) v0(end-1:-2:1,:)-1i*v0(end:-2:2,:)]];
  
  k = (2.*pi./d).*[0:Nh-1 0 -Nh+1:-1]';   % wave numbers
  L = k.^2 - k.^4;                        % Fourier multipliers
  E = exp(h*L);  E2 = exp(h*L/2);
  M = 16;                                 % no. of points for complex means
  r = exp(1i*pi*((1:M)-.5)/M);            % roots of unity
  LR = h*L(:,ones(M,1)) + r(ones(N,1),:);
  Q = h*real(mean((exp(LR/2)-1)./LR ,2));
  f1 = h*real(mean((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3 ,2));
  f2 = h*real(mean((2+LR+exp(LR).*(-2+LR))./LR.^3 ,2));
  f3 = h*real(mean((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3 ,2));
  at = a0;  tt = 0;  g = 0.5i*k*N;
  
  E =  repmat(E,1,NN+1);  E2 = repmat(E2,1,NN+1);  Q = repmat(Q,1,NN+1);
  f1 = repmat(f1,1,NN+1); f2 = repmat(f2,1,NN+1);  f3 = repmat(f3,1,NN+1);
  g = [g repmat(2.*g,1,NN)];
  for n = 1:nstp,
    t = n*h;                  rfv = real(ifft(v)); Nv = g.*fft(repmat(rfv(:,1),1,NN+1).*rfv);
    a = E2.*v + Q.*Nv;        rfv = real(ifft(a)); Na = g.*fft(repmat(rfv(:,1),1,NN+1).*rfv);
    b = E2.*v + Q.*Na;        rfv = real(ifft(b)); Nb = g.*fft(repmat(rfv(:,1),1,NN+1).*rfv);
    c = E2.*a + Q.*(2*Nb-Nv); rfv = real(ifft(c)); Nc = g.*fft(repmat(rfv(:,1),1,NN+1).*rfv);
    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3; end
  vt = zeros(N-2,NN);  vt(1:2:end,:) = real(v(2:Nh,2:end));
  vt(2:2:end,:) = imag(v(2:Nh,2:end));
  at(1:2:end) = real(v(2:Nh,1)); at(2:2:end) = imag(v(2:Nh,1));
