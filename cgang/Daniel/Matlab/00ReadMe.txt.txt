siminos/cgang/Daniel/Matlab 00ReadMe.txt
$Author: dborrero $
$Date: 2012-04-26 16:34:32 -0400 (Thurs, 26 Apr 2012) $
------------------------------------------------------------

This folder contains Matlab code pertaining to the 2-mode
SO(2) symmetric system that Chaos Gang is adapting from
Dangelmayr(1986) and  Porter and Knobloch (2005).

List of files:

EOMComplexTo4.m by Daniel Borrero 4/26/2012
------------------------------------------------
Calculates the equations of motion for the 2-mode (Dangelmayr)
system as system of 4 real valued equation in Cartesian coordinates,
a system of 2 polar coordinate pairs, and a 3D system with 2 radial
and 1 angular coordinates.

TwoModeSymmetryCheck.m by Daniel Borrero 4/26/2012
------------------------------------------------
Checks that the 2-mode equations are SO(2) equivariant
with z1 equivariant for m = 1 and z2 equivariant for the m = 2 mode
by analytically verifying that its Lie derivative vanishes everywhere.

TwoModeOriginStability.m by Daniel Borrero 4/27/2012
---------------------------------------------------
This calculates the fixed points of the 2-mode system
in 4D Cartesian coordinates. Only the origin is a fixed
point. It then calculates the stability of the origin
and shows that its stability is strictly determined by 
the signs of mu1 and mu2.

PKCart2Polar by Daniel Borrero 4/27/2012
-----------------------------------------
Converts points in the Porter-Knobloch system
from a 4D Cartesian representation (x1,y1,x2,y2) to a 4D
(2 radial + 2 angular) polar representation (r1,theta1,r2,theta2)
or a 3D (2 radial + 1 relative angular) polar representation
(r1,r2,Psi);

PKPolar2Cart by Daniel Borrero 4/27/2012
-----------------------------------------
Converts points in the Porter-Knobloch system
from a 3D (2 radial + 1 relative angular) polar representation
(r1,r2,Psi) to a 4D Cartesian representation (x1,y1,x2,y2)

PKPolarPtStability.m by Daniel Borrero 4/27/2012
-------------------------------------------------
Calculates the eigenvalues and eigenvectors of 
the stability matrix for the Porter-Knobloch system in
for a generic point given in the 4D Cartesian representation
(X1,Y1,X2,Y2) for a given array of parameters
params = {a1,a2,b1,b2,c1,c2,mu1,mu2,2}.