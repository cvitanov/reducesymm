function [ J ] = Jacobian( x0, t )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    r1 = 28;
    r2 = 0;
    b = 8/3;
    e = 1/10;
    s = 10;
    [Tout, Yout] = ode45(@ComplexLorenz, [0 t], x0,odeset('RelTol', 10^-6));
    l = length(Tout);
    I = [[1 0 0 0 0];[0 1 0 0 0];[0 0 1 0 0];[0 0 0 1 0];[0 0 0 0 1]];
    J=I;
    for i = 2:l
        x = Yout(i,:);
        A = [[-s 0 s 0 0];[0 -s 0 s 0];[r1-x(5) -r2 -1 -e -x(1)];[r2 r1-x(5) e -1 -x(2)];[x(3) x(4) x(1) x(2) -b]];
        J = J*(I + (Tout(i)-Tout(i-1))*A);
    end
end

