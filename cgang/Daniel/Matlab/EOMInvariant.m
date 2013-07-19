function dx = EOMInvariant(t,x)
dx = zeros(4,1);

mu1 = -0.38;
mu2 = 0.38;
c1 = -1;
c2 = 1;
a1 = -1.31;
a2 = -2.6;
b1 = 1.504;
b2 = 0.22;
e2 = 1.61;

u = x(1);
v = x(2);
w = x(3);
q = x(4);

% dx = [udot vdot wdot qdot]
dx(1) = c1*w + 2*u*(mu1 + a1*u + b1*v);
dx(2) = c2*w + 2*mu2*v + 2*b2*v^2 + 2*a2*u*v;
dx(3) = 2*mu1*w - e2*q + mu2*w + 2*c2*u^2 + 2*a1*u*w + a2*u*w + 4*c1*u*v + 2*b1*v*w + b2*v*w;
dx(4) = e2*w + 2*mu1*q + mu2*q + 2*b1*q*v + b2*q*v + q*u*(2*a1 + a2);