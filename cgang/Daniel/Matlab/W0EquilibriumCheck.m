% W0EquilibriumCheck.m by Daniel Borrero 5/8/2012
% Checks that the equilibrium reported by Bryce Robbins
% in the two mode blog on 4/28 is indeed a fixed point of the flow
% (after correcting some factors of two that were missing in the
% original equations of motion that Bryce used)
% -----------------------------------------------------------------

% Define symbolic variables
clc
clear all
syms a1 a2 b1 b2 c1 c2 mu1 mu2 e2 real
syms u v positive
syms w q real

% Two-mode equations of motion as derived by PKComplexToInvariant.m
% and recorded in the two mode blog
udot = c1*w + 2*u*(mu1 + a1*u + b1*v);
vdot = c2*w + 2*mu2*v + 2*b2*v^2 + 2*a2*u*v;
wdot = 2*mu1*w - e2*q + mu2*w + 2*c2*u^2 + 2*a1*u*w + a2*u*w + 4*c1*u*v + 2*b1*v*w + b2*v*w;
qdot = e2*w + 2*mu1*q + mu2*q + 2*b1*q*v + b2*q*v + q*u*(2*a1 + a2);

% Fixed point reported by Bryce on 4/28, corrected by DB on 5/8/12 for
% pesky factors of two
U = (b1*mu2 - b2*mu1)/(a1*b2 - a2*b1);
V = -(a1*mu2 - a2*mu1)/(a1*b2 - a2*b1);
W = 0;
Q = simple(2*U/e2*(2*c1*V+c2*U));

% Check that the point (U,V,W,Q) is actually a fixed point
XDOT = simple(subs([udot,vdot,wdot,qdot],{u,v,w,q},{U,V,W,Q}))
