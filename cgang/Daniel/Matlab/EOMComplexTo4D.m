% EOMComplexTo4.m by Daniel Borrero 4/26/2012
% ------------------------------------------------
% Calculates the equations of motion for the 2-mode (Dangelmayr)
% system as system of 4 real valued equation in Cartesian coordinates,
% a system of 2 polar coordinate pairs, and a 3D system with 2 radial
% and 1 angular coordinates.

clear all

% Intialize symbolic variables
syms x1 x2 y1 y2 a1 a2 b1 b2 c1 c2 mu1 mu2 e2 real
syms r1dot theta1 theta1dot r2dot theta2 theta2dot real
syms r1 r2 positive
syms Psi real

% Express complex variables z in terms of real and imaginary parts
% Notice that 1's form one complex number and 2's form another, whereas
% in the Complex Lorenz example of Chaosbook, x's form one complex number
% and y's form another.
z1 = x1 + i*y1;
z2 = x2 + i*y2;

% Define zdot's as per Siminos's eq:AGpolar in Atlas blog.
z1dot = mu1*z1 + a1*z1*z1*conj(z1) + b1*z1*z2*conj(z2) + c1*conj(z1)*z2;
z2dot = (mu2 -i*e2)*z2 + a2*z2*z1*conj(z1) + b2*z2*z2*conj(z2) + c2*z1*z1; 

% Calculate 4D flow by taking real and imaginary parts of zdot equations
disp('The 2-mode system expressed in terms of 4 real variables')
x1dot = simplify(real(z1dot))
y1dot = simplify(imag(z1dot))
x2dot = simplify(real(z2dot))
y2dot = simplify(imag(z2dot))

% Solve for equations in polar coordinates 
% {x1,y1,x2,y2} -> {r1*cos(theta1),r1*sin(theta1),r2*cos(theta2),r2*sin(theta2)}
eq1 = subs(x1dot,{x1,y1,x2,y2},{r1*cos(theta1),r1*sin(theta1),r2*cos(theta2),r2*sin(theta2)})...
      -r1dot*cos(theta1) + r1*sin(theta1)*theta1dot;
eq2 = subs(y1dot,{x1,y1,x2,y2},{r1*cos(theta1),r1*sin(theta1),r2*cos(theta2),r2*sin(theta2)})...
      -r1dot*sin(theta1) - r1*cos(theta1)*theta1dot;
eq3 = subs(x2dot,{x1,y1,x2,y2},{r1*cos(theta1),r1*sin(theta1),r2*cos(theta2),r2*sin(theta2)})...
      -r2dot*cos(theta2) + r2*sin(theta2)*theta2dot;
eq4 = subs(y2dot,{x1,y1,x2,y2},{r1*cos(theta1),r1*sin(theta1),r2*cos(theta2),r2*sin(theta2)})...
      -r2dot*sin(theta2) - r2*cos(theta2)*theta2dot;

S = solve(eq1,eq2,eq3,eq4,'r1dot','r2dot','theta1dot','theta2dot');

disp('Substituting {r1*cos(theta1),r1*sin(theta1),r2*cos(theta2),r2*sin(theta2)}')
disp('for {x1,y1,x2,y2} and solving for the velocities in these coordinates, we get:')
r1dot = simple(S.r1dot)
theta1dot = simple(S.theta1dot)
r2dot = simple(S.r2dot)
theta2dot = simple(S.theta2dot)

disp('Notice that theta1 and theta2 only appear in the combination (2*theta1 - theta2)')
disp('Recasting the equations in terms of Psi = 2*theta1 - theta2,')
disp('we can reduce the problem to a 3D system')
r1dot = subs(r1dot,2*theta1-theta2,Psi)
r2dot = subs(r2dot,2*theta1-theta2,Psi)
psidot = subs(2*theta1dot - theta2dot,2*theta1-theta2,Psi)

