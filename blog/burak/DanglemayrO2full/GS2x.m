function v = GS2x(vGS)

%Cartesian basis for full state space:

ex1 = [1;0;0;0];
ex2 = [0;1;0;0];
ey1 = [0;0;1;0];
ey2 = [0;0;0;1];

%Gram-Schmidt basis:

ex1GS = LieElement(pi/4,ex1);
ex2GS = LieElement(pi/4,ex2);
ey1GS = LieElement(pi/4,ey1);
ey2GS = LieElement(pi/4,ey2);

V = [ex1GS ex2GS ey1GS ey2GS];

vGS = [0;
	   vGS];

v = V * vGS;
