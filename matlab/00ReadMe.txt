siminos/matlab/00ReadMe.txt
$Author$
$Date:2010-04-21 18:36:47 -0400 (Wed, 21 Apr 2010) $
---------------------------------------------------------------

Ruslan 2011-03-29
My KS related Matlab code for student's enjoyment (or pain)...
	matlab/ruslan/

[BELOW: OLD STUFF, REMOVE]
ComplexLorenz.m is just the function for the complex Lorenz
equations used for ode45 in PlotCL.m and FiniteTimeStep.m.

PlotCL.m takes in the initial value and length of time you want
to integrate over, and plots the trajectory of the system.

FiniteTimeStep.m takes in the initial value, length of the time
interval before it is rotated back into the x1=0 slice, and the
number of intervals this is to be done one. It outputs the
values of the system after each rotation and displays these in
a graph (this is not a graph of the entire trajectory).

CL4D.m is the function used in ode45 for the reduced complex
Lorenz equations for Infinitesimal.m.

Infinitesimal.m takes in the initial value and length of time
you want to integrate over. It rotates the initial point into
the x1=0 hyperplane and then calculates the trajectory using
the reduced equations you get in the method of moving frames.

HilbEq.m is the function used in ode45 for a Hilbert basis
version of the system.

Hilbert1.m calculates the trajectory of the flow in the
coordinates of a Hilbert basis.

Hilbert2.m calculates the trajectory of the flow in standard
coordinates then converts the resulting flow into the
coordinates of a Hilbert basis.

---------------------------------------------------------------
Predrag 2010-04-21
was siminos/chao/matlab/00ReadMe.txt, but Tortoise could not
read the files, so Evangelos copied them to here, and svn removed
    siminos/chao/matlab/
