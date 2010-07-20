function [ T,theta,x ] = rpoNewton( T0, theta0, x0 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    x = transpose(x0);
    T = T0;
    theta = -theta0;
    I = [[1 0 0 0 0];[0 1 0 0 0];[0 0 1 0 0];[0 0 0 1 0];[0 0 0 0 1]];
    t = [[0 -1 0 0 0];[1 0 0 0 0];[0 0 0 -1 0];[0 0 1 0 0];[0 0 0 0 0]];
    d = [1];
    hold on;
    while norm(d) > 10^-12
        [f J] = Jacobian(transpose(x),T);
        v = ComplexLorenz(0,x);
        Tx = t*x;
        Tfx = t*transpose(f);
        G = [[cos(theta) sin(theta) 0 0 0];[-sin(theta) cos(theta) 0 0 0];
            [0 0 cos(theta) sin(theta) 0];[0 0 -sin(theta) cos(theta) 0];
            [0 0 0 0 1]];
        M = [[(G*J-I) v Tfx];[transpose(v) 0 0];[transpose(Tx) 0 0]];
        b= [(x-G*transpose(f));0;0];
        d = M\b;
        dT = d(6);
        dtheta = d(7);
        dx = d(1:5);
        plot([T T+dT],[theta theta+dtheta])
        T = T + dT;
        x = x + dx;
        theta = theta + dtheta;
    end
    theta = -theta;
end

