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
