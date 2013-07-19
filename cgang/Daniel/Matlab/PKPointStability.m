% PKPolarPtStability.m by Daniel Borrero 4/27/2012
% -------------------------------------------------
% Calculates the eigenvalues and eigenvectors of 
% the stability matrix for the Porter-Knobloch system in
% for a generic point given in the 4D Cartesian representation
% (X1,x2,y1,Y2) for a given array of parameters
% params = {a1,a2,b1,b2,c1,c2,mu1,mu2,2}.

function [E,Lambda] = PKPointStability(X1,x2,y1,Y2,params)
% Intialize symbolic variables
syms x1 y1 x2 y2 a1 a2 b1 b2 c1 c2 mu1 mu2 e2 real

% Specify 2-mode flow in Cartesian coordinates as solved for by EOMComplexTo4D.m
x1dot = a1*x1^3 + b1*x1*y1^2 + c1*x1*y1 + a1*x1*x2^2 + b1*x1*y2^2 + mu1*x1 + c1*x2*y2;
x2dot = a1*x1^2*x2 + c1*x1*y2 + b1*y1^2*x2 - c1*y1*x2 + a1*x2^3 + b1*x2*y2^2 + mu1*x2;
y1dot = a2*x1^2*y1 + c2*x1^2 + b2*y1^3 + a2*y1*x2^2 + b2*y1*y2^2 + mu2*y1 - c2*x2^2 + e2*y2;
y2dot = a2*x1^2*y2 + 2*c2*x1*x2 + b2*y1^2*y2 - e2*y1 + a2*x2^2*y2 + b2*y2^3 + mu2*y2;

% Define vectors of coordinates and time derivatives
X = [x1 x2 y1 y2]';
V = [x1dot x2dot y1dot y2dot]';

% Calculate the stability matrix
A = jacobian(V,X)

% Load system parameters
A1 = params{1}; A2 = params{2}; B1 = params{3}; B2 = params{4};
C1 = params{5}; C2 = params{6}; MU1 = params{7}; MU2 = params{8}; 
E2 = params{9};

% Evaluate stability matrix at point
A = subs(A,{x1,x2,y1,y2,a1,a2,b1,b2,c1,c2,mu1,mu2,e2},...
    {X1,x2,y1,Y2,A1,A2,B1,B2,C1,C2,MU1,MU2,E2});

% Calculate eigenvectors E and eigenvalues Lambda
[E,Lambda] = eig(A)

end

