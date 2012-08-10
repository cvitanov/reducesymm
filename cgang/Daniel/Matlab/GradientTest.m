clear all

syms a1 a2 b1 b2 c1 c2 mu1 mu2 e2 real
syms u v positive
syms w q real

X = [u v w q];

udot = c1*w + 2*u*(mu1 + a1*u + b1*v);
vdot = c2*w + 2*mu2*v + 2*b2*v^2 + 2*a2*u*v;
wdot = 2*mu1*w - e2*q + mu2*w + 2*c2*u^2 + 2*a1*u*w + a2*u*w + 4*c1*u*v + 2*b1*v*w + b2*v*w;
qdot = e2*w + 2*mu1*q + mu2*q + 2*b1*q*v + b2*q*v + q*u*(2*a1 + a2);

V = [udot vdot wdot qdot];

s = w^2 + q^2 - 4*u^2*v;

f = gradient(s,X)

A = jacobian(V,X)

A(1,3)