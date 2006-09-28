function wdot = ksequil(t, w)
%  Equation for KS equilibrium: u_xxxx + u_xx + u*u_x = 0
%  w = (u, u_x, u_xx, u_xxx)

  wdot = [w(2:4); -w(3)-w(2).*w(1)];
  