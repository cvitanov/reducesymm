function [ du ] = HilbEq( t, u )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    r1 = 28;
    r2 = 0;
    b = 8/3;
    e = 1/10;
    sig = 10;
    du(1,1) = 2*sig*(u(4)-u(1));
    du(2,1) = -2*u(2)-2*u(4)*(u(5)-r1)+2*r2*u(3);
    du(3,1) = e*u(4)-(sig+1)*u(3)+r2*u(1);
    du(4,1) = sig*u(2)-(sig+1)*u(4)-e*u(3)+u(1)*(r1-u(5));
    du(5,1) = u(4)-b*u(5);

end

