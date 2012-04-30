% PKComplexToInvariant.m by Daniel Borrero 4/29/2012
% ------------------------------------------------
% This code verifies the various transformations
% between the complex formulation of the Porter-Knobloch
% system and its formulation in terms of the invariant
% bases u, v, w, and q as described in section s:twoMode
% of the 2mode project/report.

clc
clear all
syms a1 a2 b1 b2 c1 c2 mu1 mu2 e2 real
syms z1 z2
syms u v positive
syms w q real
syms theta1 theta2 real
syms x1 y1 x2 y2 real
syms phi

%Define z's in terms of u,v, and theta's
z1 = sqrt(u)*(cos(theta1) + i*sin(theta1)); 
z2 = sqrt(v)*(cos(theta2) + i*sin(theta2));

W = simple(z1^2*conj(z2) + conj(z1)^2*z2);
W2 = 2*real(z1^2*conj(z2)); % w in terms of real parts of z1^2*conj(z2)
DW = simple(W-W2) %Check Dang86(1.2)polar is correct.
W = subs(W,2*theta1-theta2,phi)%Check Dang86(1.2)polar is correct.

Q = simple(z1^2*conj(z2) - conj(z1)^2*z2)/i;
Q2 = 2*imag(z1^2*conj(z2));  % q in terms of imaginary parts of z1^2*conj(z2)
DQ = simple(Q-Q2)%Check Dang86(1.2)polar is correct.
Q = subs(Q,2*theta1-theta2,phi)%Check Dang86(1.2)polar is correct.

C = simple((W/2)^2+(Q/2)^2-u^2*v) % Check syzygy

%Define zdot's
z1dot = mu1*z1 + a1*z1*z1*conj(z1) + b1*z1*z2*conj(z2) + c1*conj(z1)*z2;
z2dot = (mu2 -i*e2)*z2 + a2*z2*z1*conj(z1) + b2*z2*z2*conj(z2) + c2*z1*z1; 

% Check expression for udot in terms of (u,v,w,q)
udot = simple(conj(z1)*z1dot + z1*conj(z1dot)); % Plug in z's into udot
udot = subs(udot,2*theta1-theta2,phi); %substitute phi = 2*theta1-theta2
udot = simple(subs(udot,{cos(phi),sin(phi)},{w/(2*u*v^(1/2)),q/(2*u*v^(1/2))})) % Substitute u,v,w,q's for sin(phi)'s and cos(phi)'s

% Check expression for vdot in terms of (u,v,w,q)
vdot = simple(conj(z2)*z2dot + z2*conj(z2dot)); % Plug in z's into vdot
vdot = subs(vdot,2*theta1-theta2,phi);  %substitute phi = 2*theta1-theta2
vdot = simple(subs(vdot,{cos(phi),sin(phi)},{w/(2*u*v^(1/2)),q/(2*u*v^(1/2))})) % Substitute u,v,w,q's for sin(phi)'s and cos(phi)'s

% Check expression for wdot in terms of (u,v,w,q)
wdot = simple(2*conj(z2)*z1*z1dot + 2*z2*conj(z1)*conj(z1dot) + z1^2*conj(z2dot) + conj(z1)^2*z2dot);% Plug in z's into wdot
wdot = subs(wdot,2*theta1-theta2,phi); %substitute phi = 2*theta1-theta2
wdot = simple(subs(wdot,{cos(phi),sin(phi)},{w/(2*u*v^(1/2)),q/(2*u*v^(1/2))})) % Substitute u,v,w,q's for sin(phi)'s and cos(phi)'s

% Check expression for qdot in terms of (u,v,w,q)
qdot = simple((2*conj(z2)*z1*z1dot - 2*z2*conj(z1)*conj(z1dot) + z1^2*conj(z2dot) - conj(z1)^2*z2dot)/i);% Plug in z's into qdot
qdot = subs(qdot,2*theta1-theta2,phi); %substitute phi = 2*theta1-theta2
qdot = simple(subs(qdot,{cos(phi),sin(phi)},{w/(2*u*v^(1/2)),q/(2*u*v^(1/2))})) % Substitute u,v,w,q's for sin(phi)'s and cos(phi)'s
