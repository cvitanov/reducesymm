function [  ] = Infinitesimal( y0, tspan )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    theta = atan2(y0(1), y0(2));
    G = [[cos(-theta),sin(-theta),0,0,0];[-sin(-theta),cos(-theta),0,0,0];
        [0,0,cos(-theta),sin(-theta),0];[0,0,-sin(-theta),cos(-theta),0];
        [0,0,0,0,1]];
    y0 = G * transpose(y0);
    y0 = transpose(y0);
    y0 = y0(2:5);
    [Tout, Yout] = ode45(@CL4D, [0 tspan], y0,odeset('RelTol', 10^-6));
    plot3(Yout(:,1), Yout(:,2), Yout(:,3))
    l = length(Yout(:,1));
    y = Yout(l,:)
end

