function [  ] = PlotCL(y0, tspan)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    [Tout, Yout] = ode45(@ComplexLorenz, [0 tspan], y0,odeset('RelTol', 10^-6));
    plot3(Yout(:,1), Yout(:,2), Yout(:,3))
end

