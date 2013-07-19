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
and 1 angular coordinates. Also calculates stability matrix and its
trace in the {x1,x2,y1,y2} basis.

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

PKComplexToInvariant.m by Daniel Borrero 4/29/2012
------------------------------------------------
This code verifies the various transformations
between the complex formulation of the Porter-Knobloch
system and its formulation in terms of the invariant
bases u, v, w, and q as described in section s:twoMode
of the 2mode project/report. Also calculates the stability
matrix and its trace in the {u,v,w,q} basis.

W0EquilibriumCheck.m by Daniel Borrero 5/8/2012
--------------------------------------------------
Checks that the equilibrium (expressed in the invariant
polynomial basis) reported by Bryce Robbins 
in the two mode blog on 4/28 is indeed a fixed point of the flow
(after correcting some factors of two that were missing in the
original equations of motion that Bryce used)

Sobolev_product.m by Daniel Borrero 08/3/2012
--------------------------------------------------
Calculates the inner product of two vectors for the 2-mode
system of Porter & Knobloch using the Sobolev H^1 norm whose 
metric tensor is given by

                   [[1 0 0 0];
               g =  [0 1 0 0];
                    [0 0 4 0];
                    [0 0 0 4]]

so that <x|x'> = x1*x1' + x2*x2' + 4(x3*x3' + x4*x4')

PKRoots.m by Daniel Borrero 08/08/2012
-----------------------------------------------------------
Finds the roots for the Porter & Knobloch system in the
(u,v,w,q) invariant polynomial basis for a given
set of parameters and displays them one by one, checking if they
satisfy the syzygy w^2+q^2-4*u^2*v = 0