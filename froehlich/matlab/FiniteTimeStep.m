function [y] = FiniteTimeStep( y0, tstep, iter )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    theta = atan2(y0(1), y0(2));
    G = [[cos(-theta),sin(-theta),0,0,0];[-sin(-theta),cos(-theta),0,0,0];
        [0,0,cos(-theta),sin(-theta),0];[0,0,-sin(-theta),cos(-theta),0];
        [0,0,0,0,1]];
    a = G * transpose(y0);
    y(1,:) = transpose(a);
    for i = 1:iter
        [Tout, Yout] = ode45(@ComplexLorenz, [0, tstep],y(i,:),odeset('RelTol', 10^-6));
        l = length(Yout(:,1));
        y(i+1,:) = Yout(l,:);
        theta = atan2(y(i+1,1),y(i+1,2));
        G = [[cos(-theta),sin(-theta),0,0,0];[-sin(-theta),cos(-theta),0,0,0];
            [0,0,cos(-theta),sin(-theta),0];[0,0,-sin(-theta),cos(-theta),0];
            [0,0,0,0,1]];
        a = G * transpose(y(i+1,:));
        y(i+1,:) = transpose(a);
    end
    plot3(y(:,1),y(:,2), y(:,3))
end

