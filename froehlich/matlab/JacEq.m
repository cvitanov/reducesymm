function [ dy ] = JacEq( t, y )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    r1 = 28;
    r2 = 0;
    b = 8/3;
    e = 1/10;
    sig = 10;
    dy(1) = (-sig*y(1))+ sig*y(3);
    dy(2) = (-sig*y(2))+ sig*y(4);
    dy(3) = (r1 - y(5))*y(1) - r2*y(2) - y(3) - e*y(4);
    dy(4) = r2*y(1) + (r1 - y(5))*y(2) + e*y(3) - y(4);
    dy(5) = -b*y(5) + y(1)*y(3) + y(2)*y(4);
    J = [transpose(y(6:10));transpose(y(11:15));transpose(y(16:20));transpose(y(21:25));transpose(y(26:30))];
    A = [[-sig 0 sig 0 0];[0 -sig 0 sig 0];[r1-y(5) -r2 -1 -e -y(1)];[r2 r1-y(5) e -1 -y(2)];[y(3) y(4) y(1) y(2) -b]];
    dJ = A*J;
    dy = [dy dJ(1,:) dJ(2,:) dJ(3,:) dJ(4,:) dJ(5,:)];
    dy = transpose(dy);
end

