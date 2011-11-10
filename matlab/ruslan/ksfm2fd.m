function [afd, ell, rho] = ksfm2fd(a, L)
% function [afd, ell, rho] = ksfm2fd(a, L)
%   Map KSE solutions to the fundamental domain

  va = a(1:2:end,:) + 1i*a(2:2:end,:);  theta1 = angle(va(1,:));  fa = ksfm2(a,L);
  thetadot1 = (a(1,:).*fa(2,:)-a(2,:).*fa(1,:))./(a(1,:).^2+a(2,:).^2);
  rho = ones(1,size(a,2)); rho(thetadot1 > 0) = -1;
  tau = exp(1i*(1:size(va,1))'*(pi/2-theta1));  vf = tau.*va;
  ell = mod((pi/2 - theta1)*(L/2/pi)+L/2,L)-L/2;
  afd = zeros(size(a));  afd(1:2:end,:) = real(vf);  afd(2:2:end,:) = imag(vf);
  ir = find(rho == -1);  afd(1:2:end,ir) = -afd(1:2:end,ir);
