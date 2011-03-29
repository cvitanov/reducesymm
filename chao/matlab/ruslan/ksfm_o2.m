function [ao, s, r] = ksfm_o2(a, L, k)
% function [as, s] = ksfm_o2(a, L, k)
%   Map KSE solutions to M/O(2) using the phase of 1st Fourier mode (k=0),
%   or the difference between phases k+1 and k (k > 0)

  if nargin < 3, k = 0; end
  va = a(1:2:end,:) + 1i*a(2:2:end,:);  theta1 = angle(va(k+1,:));
  fa = ksfm2(a, L);  k1 = 2*k+1; k2 = 2*k+2;
  thetadot1 = (a(k1,:).*fa(k2,:)-a(k2,:).*fa(k1,:))./(a(k1,:).^2+a(k2,:).^2);
  if k == 0, theta2 = pi/2; thetadot2 = 0; 
  else theta2 = angle(va(k,:)); k1 = 2*k-1; k2 = 2*k;
    thetadot2 = (a(k1,:).*fa(k2,:)-a(k2,:).*fa(k1,:))./(a(k1,:).^2+a(k2,:).^2); end
  s = mod((theta2 - theta1)/2/pi+0.5, 1) - 0.5;
  tau = exp(2i*pi*(1:size(va,1))'*s);  vf = tau.*va;
  ao = zeros(size(a));  ao(1:2:end,:) = real(vf);  ao(2:2:end,:) = imag(vf);
  r = ones(1,size(a,2)); r(thetadot1-thetadot2 > 0) = -1;
  ir = find(r == -1);  ao(1:2:end,ir) = -ao(1:2:end,ir);
