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