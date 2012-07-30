% TwoModeOriginStability.m by Daniel Borrero 4/27/2012
% ---------------------------------------------------
% This calculates the fixed points of the 2-mode system
% in 4D Cartesian coordinates. Only the origin is a fixed
% point. It then calculates the stability of the origin
% and shows that its stability is strictly determined by 
% the signs of mu1 and mu2.

clear all

% Intialize symbolic variables
syms x1 y1 x2 y2 a1 a2 b1 b2 c1 c2 mu1 mu2 e2 real

% Specify 2-mode flow in Cartesian coordinates as solved for by EOMComplexTo4D.m
x1dot = a1*x1^3 + b1*x1*y1^2 + c1*x1*y1 + a1*x1*x2^2 + b1*x1*y2^2 + mu1*x1 + c1*x2*y2;
x2dot = a1*x1^2*x2 + c1*x1*y2 + b1*y1^2*x2 - c1*y1*x2 + a1*x2^3 + b1*x2*y2^2 + mu1*x2;
y1dot = a2*x1^2*y1 + c2*x1^2 + b2*y1^3 + a2*y1*x2^2 + b2*y1*y2^2 + mu2*y1 - c2*x2^2 + e2*y2;
y2dot = a2*x1^2*y2 + 2*c2*x1*x2 + b2*y1^2*y2 - e2*y1 + a2*x2^2*y2 + b2*y2^3 + mu2*y2;

S1 = solve(x1dot,x2dot,y1dot,y2dot,'x1','x2','y1','y2');
disp('The fixed points of the CL system are:')
EQ = double([S1.x1, S1.x2, S1.y1, S1.y2]) % This only returns the origin

% Define vectors of coordinates and time derivatives
disp('First, we define the vector of coordinates X and the flow V.')
X = [x1 x2 y1 y2]'
V = [x1dot x2dot y1dot y2dot]';

% Calculate the stability matrix
disp('Then, we calculate stability matrix, A = jacobian(V,X) at {0,0,0,0}')
A = jacobian(V,X)
%tr = simple(trace(A) - (2*(mu1 + mu2) + (4*a1+2*a2)*(x1^2+x2^2) + (2*b1+ 4*b2)*(y1^2+y2^2)))
tr = simple(trace(A))

A = subs(A,{x1,x2,y1,y2},{0,0,0,0})

disp('The eigenvalues and eigenvectors of the stability matrix at the origin are')

[E,Lambda] = eig(A)