function [ dy ] = CL4D( t, y )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    r1 = 28;
    r2 = 0;
    b = 8/3;
    e = 1/10;
    sig = 10;
    dy(1,1) = sig*(y(3) - y(1));
    dy(2,1) = -y(2)+r2*y(1)-(e + sig*y(2)/y(1))*y(3);
    dy(3,1) = -y(3) + (r1 - y(4))*y(1) + (e + sig*y(2)/y(1))*y(2);
    dy(4,1) = -b*y(4) + y(1)*y(3);

end

