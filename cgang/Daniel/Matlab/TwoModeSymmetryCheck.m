% TwoModeSymmetryCheck.m by Daniel Borrero 4/26/2012
% ------------------------------------------------
% Checks that the 2-mode equations are SO(2) equivariant
% with z1 equivariant for m = 1 and z2 equivariant for the m = 2 mode
% by analytically verifying that its Lie derivative vanishes everywhere.

function TwoModeSymmetryCheck
clear all

% Intialize symbolic variables
syms x1 y1 x2 y2 a1 a2 b1 b2 c1 c2 mu1 mu2 e2 real

% Specify 2-mode flow in Cartesian coordinates as solved for by EOMComplexTo4D.m
x1dot = a1*x1^3 + b1*x1*y1^2 + c1*x1*y1 + a1*x1*x2^2 + b1*x1*y2^2 + mu1*x1 + c1*x2*y2;
x2dot = a1*x1^2*x2 + c1*x1*y2 + b1*y1^2*x2 - c1*y1*x2 + a1*x2^3 + b1*x2*y2^2 + mu1*x2;
y1dot = a2*x1^2*y1 + c2*x1^2 + b2*y1^3 + a2*y1*x2^2 + b2*y1*y2^2 + mu2*y1 - c2*x2^2 + e2*y2;
y2dot = a2*x1^2*y2 + 2*c2*x1*x2 + b2*y1^2*y2 - e2*y1 + a2*x2^2*y2 + b2*y2^3 + mu2*y2;

% Define vectors of coordinates and time derivatives
disp('First, we define the vector of coordinates X and the flow V.')
X = [x1 x2 y1 y2]'
V = [x1dot x2dot y1dot y2dot]';

% Calculate the stability matrix
disp('Then, we calculate stability matrix, A = jacobian(V,X)')
A = jacobian(V,X)

% Create the 2-mode Lie group generator
disp('The Lie group generator T for the 2-mode system is given by')
T = zeros(4,4);
T(2,1) = -1;
T(1,2) = 1;
T(4,3) = -2;
T(3,4) = 2

% Check that the Lie derivative vanishes
disp('Now check that Lie derivative (T*v - A*T*X) vanishes')
simplify(T*V - A*T*X)

disp('This vanishes, so this system of equations has SO(2) symmetry with 2 modes')
end