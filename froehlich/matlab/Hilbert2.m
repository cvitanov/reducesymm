function [  ] = Hilbert2( y0, tspan )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [Tout, Yout] = ode45(@ComplexLorenz, [0 tspan], y0,odeset('RelTol', 10^-6));
    l = length(Yout(:,1));
    for i = 1:l
        u(i,1) = Yout(i,1)^2+Yout(i,2)^2;
        u(i,2) = Yout(i,3)^2+Yout(i,4)^2;
        u(i,3) = Yout(i,1)*Yout(i,4)-Yout(i,2)*Yout(i,3);
        u(i,4) = Yout(i,1)*Yout(i,3)+Yout(i,2)*Yout(i,4);
        u(i,5) = Yout(i,5);
    end
    plot3(u(:,1), u(:,2), u(:,3))

end

