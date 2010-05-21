function [ dy ] = ComplexLorenz( t, y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    r1 = 28;
    r2 = 0;
    b = 8/3;
    e = 1/10;
    sig = 10;
    dy(1,1) = (-sig*y(1))+ sig*y(3);
    dy(2,1) = (-sig*y(2))+ sig*y(4);
    dy(3,1) = (r1 - y(5))*y(1) - r2*y(2) - y(3) - e*y(4);
    dy(4,1) = r2*y(1) + (r1 - y(5))*y(2) + e*y(3) - y(4);
    dy(5,1) = -b*y(5) + y(1)*y(3) + y(2)*y(4);

end

