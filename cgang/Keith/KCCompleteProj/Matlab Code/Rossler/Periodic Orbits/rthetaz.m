% Written by K. Carroll 5/1/2012
function [ k ] = rthetaz( xpos )
    a = 0.2;
    b = 0.2;
    c = 5.7;
    r = xpos(1:1,1);
    theta = xpos(1:1,2);
    z = xpos(1:1,3);
    t2 = 1/(1+z/r*sin(theta)+a/2*sin(2*theta));
    
    k(1:1,1) = (-z*cos(theta)+a*r*(sin(theta))^2)*t2;
    k(1:1,2) = 1;
    k(1:1,3) = (b+z*(r*cos(theta)-c))*t2;    
    k(1:1,4) = t2;
end