function [ x,J ] = Jacobian3( x0, t )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    r1 = 28;
    r2 = 0;
    b = 8/3;
    e = 1/10;
    s = 10;
    J0 = [1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1];
    y0 = [x0 J0];
    [Tout, Yout] = ode45(@JacEq, [0 t], y0,odeset('RelTol', 10^-12));
    l = length(Tout);
    x = Yout(l,1:5);
    J = [Yout(l,6:10);Yout(l,11:15);Yout(l,16:20);Yout(l,21:25);Yout(l,26:30)];
end

