% This function defines the extended ODE system for the 
% Porter Knobloch system in invariant polynomial basis (u,v,w,q),
% for calculation of its Lyapunov exponents using the method of Wolf 
% et al. Physica D Vol. 16, pp. 285-317, 1985
function OUT = PKeqInv(t,X,params)

% Define parameters
mu1 = params(1);
mu2 = params(2);
a1 = params(3);
a2 = params(4);
b1 = params(5);
b2 = params(6);
c1 = params(7);
c2 = params(8);
e2 = params(9);

% Define coordinates
u = X(1);
v = X(2);
w = X(3);
q = X(4);

% Evaluate time derivatives
udot = c1*w + 2*u*(mu1 + a1*u + b1*v); % udot
vdot = c2*w + 2*mu2*v + 2*b2*v^2 + 2*a2*u*v; % vdot
wdot = 2*mu1*w - e2*q + mu2*w + 2*c2*u^2 + 2*a1*u*w + a2*u*w + 4*c1*u*v + 2*b1*v*w + b2*v*w; % wdot
qdot = e2*w + 2*mu1*q + mu2*q + 2*b1*q*v + b2*q*v + q*u*(2*a1 + a2); % qdot

% Q is used for variational solver and is set by lyapunov.m. Ignore.
Q = [[X(5), X(9),  X(13), X(17)];
     [X(6), X(10), X(14), X(18)];
     [X(7), X(11), X(15), X(19)];
     [X(8), X(12), X(16), X(20)]];

% Evaluate Jacobian
J = [[2*mu1 + 4*a1*u + 2*b1*v, 2*b1*u, c1, 0];...
[2*a2*v, 2*mu2 + 2*a2*u + 4*b2*v, c2, 0];...
[2*a1*w + a2*w + 4*c2*u + 4*c1*v, 4*c1*u + 2*b1*w + b2*w, 2*mu1 + mu2 + 2*a1*u + a2*u + 2*b1*v + b2*v, -e2];...
[q*(2*a1 + a2), 2*b1*q + b2*q, e2, 2*mu1 + mu2 + 2*b1*v + b2*v + u*(2*a1 + a2)]];

% Set up variational equation
F = J*Q;

% Package output for use with lyapunov.m
OUT = [udot; vdot; wdot; qdot; F(:)];