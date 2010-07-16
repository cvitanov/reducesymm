function [ ] = ReducedAttractor( x,t )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [Tout, Yout] = ode45(@ComplexLorenz, [0 t], x,odeset('RelTol', 10^-6));
    l = length(Tout);
    y0 = Yout(l,:);
    [T1, Y1] = ode45(@ComplexLorenz, [0 t], y0,odeset('RelTol', 10^-6));
    T = [[0 -1 0 0 0];[1 0 0 0 0];[0 0 0 -1 0];[0 0 1 0 0];[0 0 0 0 0]];
    l=length(T1);
    Y = [];
    for i=1:l
        y = transpose(Y1(i,:));
        theta = atan2(y(1),y(2));
        G = [[cos(-theta),sin(-theta),0,0,0];[-sin(-theta),cos(-theta),0,0,0];
            [0,0,cos(-theta),sin(-theta),0];[0,0,-sin(-theta),cos(-theta),0];
            [0,0,0,0,1]];
        y = G*y;
        Y(:,i)=y;
    end
    plot3(Y(2,:),Y(4,:),Y(5,:))
end

