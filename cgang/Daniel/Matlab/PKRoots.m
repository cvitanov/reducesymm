% PKRoots.m by Daniel Borrero 08/08/2012
% -----------------------------------------------------------
% Finds the roots for the Porter & Knobloch system in the
% (u,v,w,q) invariant polynomial basis for a given
% set of parameters and displays them one by one, checking if they
% satisfy the syzygy w^2+q^2-4*u^2*v = 0

clear all
syms u v w q

% Define parameters
mu1 = -0.38;
mu2 = 0.38;
c1 = -1;
c2 = 1;
a1 = -1.31;
a2 = -2.6;
b1 = 1.504;
b2 = 0.22;
e2 = 1.61;

% Siminos 08/10/2012
mu1 = -0.234;
mu2 = 0.28;
c1 = 1;
c2 = 1;
a1 = -1.4;
a2 = -2.18;
b1 = -0.1;
b2 = 0.22;
e2 = 1.217;

eq1 = c1*w + 2*u*(mu1 + a1*u + b1*v); % udot
eq2 = c2*w + 2*mu2*v + 2*b2*v^2 + 2*a2*u*v; % vdot
eq3 = 2*mu1*w - e2*q + mu2*w + 2*c2*u^2 + 2*a1*u*w + a2*u*w + 4*c1*u*v + 2*b1*v*w + b2*v*w; % wdot
eq4 = e2*w + 2*mu1*q + mu2*q + 2*b1*q*v + b2*q*v + q*u*(2*a1 + a2); % qdot

% Solve (udot,vdot,wdot,qdot) = 0%
X = solve(eq1,eq2,eq3,eq4);

% Display roots
disp('The roots are:')
for i = 1:length(X.u)
    root = double([X.u(i), X.v(i), X.w(i), X.q(i)])
    syzygy = double(X.w(i)^2 + X.q(i)^2 -4*X.u(i)^2*X.v(i))
end