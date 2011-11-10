function [ar, r, irp, rrp, irm, rrm] = ksfm_rpm(a, L, k)
% function [ar, r] = ksfm2fd(a, L, k)
%   Map KSE solutions to fundamental domain w.r.t. reflection symmetry.
%     ar - reflected state, r = \pm 1 (reflection parameter), 
%     ir - index of r change, rr = (tr - tt(ir))/h - reflection time

  if nargin < 3, k = 0; end,  fa = ksfm2(a, L);  k1 = 2*k+1; k2 = 2*k+2;
  thetadot1 = (a(k1,:).*fa(k2,:)-a(k2,:).*fa(k1,:))./(a(k1,:).^2+a(k2,:).^2);
  if k == 0, thetadot2 = 0;  else k1 = 2*k-1; k2 = 2*k;
    thetadot2 = (a(k1,:).*fa(k2,:)-a(k2,:).*fa(k1,:))./(a(k1,:).^2+a(k2,:).^2); end
  dthdot = thetadot1-thetadot2;  r = ones(1,size(a,2))/2; r(dthdot > 0) = -.5;
  irp = find((dthdot(1:end-1) <= 0) & (dthdot(2:end) > 0)); 
  irm = find((dthdot(1:end-1) > 0) & (dthdot(2:end) <= 0)); 
  rrp = -dthdot(irp)./(dthdot(irp+1)-dthdot(irp));
  rrm = -dthdot(irm)./(dthdot(irm+1)-dthdot(irm));
  irr = find(r == -1);  ar = a;  ar(1:2:end,irr) = -ar(1:2:end,irr);
