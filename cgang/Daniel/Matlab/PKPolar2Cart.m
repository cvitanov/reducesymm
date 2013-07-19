% PKPolar2Cart by Daniel Borrero 4/27/2012
% -----------------------------------------
% Converts points in the Porter-Knobloch system
% from a 3D (2 radial + 1 relative angular) polar representation
% (r1,r2,Psi) to a 4D Cartesian representation (x1,x2,y1,y2)

function [x1,x2,y1,y2] = PKPolar2Cart(r1,r2,Psi)
theta2 = 0;
theta1 = (Psi + theta2)/2;
x1 = r1*cos(theta1);
x2 = r1*sin(theta1);
y1 = r2*cos(theta2);
y2 = r2*sin(theta2);
end