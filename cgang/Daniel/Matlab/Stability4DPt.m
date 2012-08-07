% 4DStability.m by Daniel Borrero 8/07/2012
% ---------------------------------------------------
% This calculates the stability of a point for the 2-mode system
% given in 4D Cartesian coordinates.

function Stability4DPt(x)

load params.mat

% Intialize symbolic variables
syms x1 y1 x2 y2 real

% Specify 2-mode flow in Cartesian coordinates as solved for by EOMComplexTo4D.m
x1dot = a1*x1^3 + a1*x1*x2^2 + b1*x1*y1^2 + c1*x1*y1 + b1*x1*y2^2 + mu1*x1 + c1*x2*y2;
x2dot = a1*x1^2*x2 + c1*x1*y2 + a1*x2^3 + b1*x2*y1^2 - c1*x2*y1 + b1*x2*y2^2 + mu1*x2;
y1dot = a2*x1^2*y1 + c2*x1^2 + a2*x2^2*y1 - c2*x2^2 + b2*y1^3 + b2*y1*y2^2 + mu2*y1 + e2*y2;
y2dot = a2*x1^2*y2 + 2*c2*x1*x2 + a2*x2^2*y2 + b2*y1^2*y2 - e2*y1 + b2*y2^3 + mu2*y2;

% Define vectors of coordinates and time derivatives
disp('First, we define the vector of coordinates X and the flow V.')
X = [x1 x2 y1 y2]'
V = [x1dot x2dot y1dot y2dot]';

% Calculate the stability matrix
disp('Then, we calculate stability matrix, A = jacobian(V,X)')
A = jacobian(V,X)

A = subs(A,{x1,x2,y1,y2},{x(1),x(2),x(3),x(4)})

disp('The eigenvalues and eigenvectors of the stability matrix are')

[E,Lambda] = eig(A)