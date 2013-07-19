% ConvInvTo4D.m 
% by Daniel Borrero 2012-08-07
% ----------------------------------------------------------------
% Converts a point given in invariant coordinates to 4D cartesian
% coordinates.

function [x1,x2,y1,y2] = ConvInvTo4D(u,v,w,q)

% First, convert point from invariant polynomial basis to polar basis
r1 = sqrt(u);
r2 = sqrt(v);
Phi = acos(w/(2*u*abs(v)^(1/2))) % Relative phase
Phi2 = asin(q/(2*u*abs(v)^(1/2)))
% Use freedom to choose overall phase and set theta2 to zero
theta1 = Phi/2;
theta2 = 0;

% Convert polar coordinates to 4D cartesian coordinates
x1 = r1*cos(theta1);
x2 = r1*sin(theta1);
y1 = r2*cos(theta2);
y2 = r2*sin(theta2);