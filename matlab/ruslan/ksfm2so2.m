function [as, s] = ksfm2so2(a, L, k)
% function [as, s] = ksfm2so2(a, L, k)
%   Map KSE solutions to M/SO(2) using the phase of 1st Fourier mode (k=0),
%   or the difference between phases k+1 and k (k > 0)

  if nargin < 3, k = 0; end
  va = a(1:2:end,:) + 1i*a(2:2:end,:);  theta1 = angle(va(k+1,:));
  if k == 0, theta2 = pi/2; else theta2 = angle(va(k,:)); end
  s = mod((theta2 - theta1)/2/pi+0.5, 1) - 0.5;
  tau = exp(2i*pi*(1:size(va,1))'*s);  vf = tau.*va;
  as = zeros(size(a));  as(1:2:end,:) = real(vf);  as(2:2:end,:) = imag(vf);
