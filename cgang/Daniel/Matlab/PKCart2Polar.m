% PKCart2Polar by Daniel Borrero 4/27/2012
% -----------------------------------------
% Converts points in the Porter-Knobloch system
% from a 4D Cartesian representation (x1,x2,y1,y2) to a 4D
% (2 radial + 2 angular) polar representation (r1,theta1,r2,theta2)
% or a 3D (2 radial + 1 relative angular) polar representation
% (r1,r2,Psi);

function [r1,theta1,r2,theta2,Psi] = PKCart2Polar(x1,x2,y1,y2)
r1 = sqrt(x1.^2+x2.^2);
r2 = sqrt(y1.^2+y2.^2);
theta1 = atan(x2/x1);
theta2 = atan(y2/y1);
Psi = 2*theta1 - theta;
end

