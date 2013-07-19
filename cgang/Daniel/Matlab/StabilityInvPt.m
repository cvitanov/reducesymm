% StabilityInvPt.m by Daniel Borrero 8/07/2012
% ---------------------------------------------------
% This calculates the stability of a point for the 2-mode system
% given in 4D Cartesian coordinates.



% Intialize symbolic variables
syms u v positive
syms w q real
syms a1 a2 b1 b2 c1 c2 mu1 mu2 e2 real

% % Define parameters
% mu1 = -0.38;
% mu2 = 0.38;
% c1 = -1;
% c2 = 1;
% a1 = -1.31;
% a2 = -2.6;
% b1 = 1.504;
% b2 = 0.22;
% e2 = 1.61;

% Specify 2-mode flow in Cartesian coordinates as solved for by EOMComplexTo4D.m
udot = c1*w + 2*u*(mu1 + a1*u + b1*v);
vdot = c2*w + 2*mu2*v + 2*b2*v^2 + 2*a2*u*v;
wdot = 2*mu1*w - e2*q + mu2*w + 2*c2*u^2 + 2*a1*u*w + a2*u*w + 4*c1*u*v + 2*b1*v*w + b2*v*w;
qdot = e2*w + 2*mu1*q + mu2*q + 2*b1*q*v + b2*q*v + q*u*(2*a1 + a2);

% Define vectors of coordinates and time derivatives
disp('First, we define the vector of coordinates X and the flow V.')
X = [u v w q];
V = [udot vdot wdot qdot]';

% Calculate the stability matrix
disp('Then, we calculate stability matrix, A = jacobian(V,X)')
A = jacobian(V,X)

%A = subs(A,{x1,x2,y1,y2},{x(1),x(2),x(3),x(4)})

disp('The eigenvalues and eigenvectors of the stability matrix are')

[E,Lambda] = eig(A)