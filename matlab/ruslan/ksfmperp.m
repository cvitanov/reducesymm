function f = ksfmperp(a, L, G)
% function f = ksfmperp(a, L, G)
%  KSFM flow normal to the SO(2) group action

  v = ksfm(0, a, L);
  t = G*a;  P = (t*t')./(t'*t);
  f = v - P*v;
  