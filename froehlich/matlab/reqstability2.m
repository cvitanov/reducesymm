function [ M ] = reqstability2( theta )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    r1 = 28;
    r2 = 0;
    b = 8/3;
    e = 1/10;
    s = 10;
    d = 1 + (e^2)/(s+1)^2;
    x = [sqrt(b*(r1-d)/d);-sqrt(b*(r1-d)*(d-1)/d);sqrt(b*d*(r1-d));0;r1-d];
    G = [[cos(theta),sin(theta),0,0,0];[-sin(theta),cos(theta),0,0,0];
        [0,0,cos(theta),sin(theta),0];[0,0,-sin(theta),cos(theta),0];
        [0,0,0,0,1]];
    x = G*x;
    A = [[-s 0 s 0 0];[0 -s 0 s 0];[r1-x(5) -r2 -1 -e -x(1)];[r2 r1-x(5) e -1 -x(2)];[x(3) x(4) x(1) x(2) -b]];
    %v = [-s*x(1)+s*x(3);-s*x(2)+s*x(4);(r1-x(5))*x(1)-r2*x(2)-x(3)-e*x(4);r2*x1+(r1-x(5))*x(2)+e*x(3)-x(4);-b*x(5)+x(1)*x(3)+x(2)*x(4)];
    T = [[0 1 0 0 0];[-1 0 0 0 0];[0 0 0 1 0];[0 0 -1 0 0];[0 0 0 0 0]];
    t = T*x;
    c = -e*s/(s+1);
    for i = 1:5
        for j=1:5
            dv = [A(1,j) A(2,j) A(3,j) A(4,j) A(5,j)];
            dT = [T(1,j) T(2,j) T(3,j) T(4,j) T(5,j)];
            Tx = T*x;
            Tx = [Tx(1) Tx(2) Tx(3) Tx(4) Tx(5)];
            M(i,j) = A(i,j)-c*T(i,j)-Tx(i)*((dv*t)/(Tx*t)-c*(dT*t)/(Tx*t));
        end
    end

end

