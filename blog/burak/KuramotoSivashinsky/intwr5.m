function [t,y] = intwr5(t, x0)

%
%Solves rosslerm using Radau 5 algorithm wrapped in odepkg
%

pkg load odepkg
          
vopt = odeset ("InitialStep", 1e-3, "MaxStep", 1e-1, "Refine", 5); %, "OutputFcn", @odephas3);
%vopt = odeset ("InitialStep", 1e-3, "MaxStep", 1e-1, "RelTol", 1e-8, "AbsTol", 1e-6, "NormControl", "on");
%vopt = odeset ("RelTol", 1e-3, "AbsTol", 1e-3, \
 %             "NormControl", "on", "OutputFcn", @odeplot);

[t,y] = ode5r (@vKS, t, x0, vopt);
%ode45 (@KSm, [0 10], randn(14,1), vopt);
