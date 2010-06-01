function [  ] = Hilbert1( y0, tspan )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    u0(1) = y0(1)^2+y0(2)^2;
    u0(2) = y0(3)^2+y0(4)^2;
    u0(3) = y0(1)*y0(4)-y0(2)*y0(3);
    u0(4) = y0(1)*y0(3)+y0(2)*y0(4);
    u0(5) = y0(5);
    [Tout Uout] = ode45(@HilbEq, [0 tspan], u0,odeset('RelTol', 10^-6));
    plot3(Uout(:,1),Uout(:,2),Uout(:,3))
    
end

