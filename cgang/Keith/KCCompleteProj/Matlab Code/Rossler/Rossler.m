%Written by K. Carroll 5/1/2012
function [k] = Rossler(xpos)
    a = 0.2;
    b = 0.2;
    c = 5.7;
    k(1:1,1) = -xpos(1:1,2)-xpos(1:1,3);
    k(1:1,2) = xpos(1:1,1)+a*xpos(1:1,2);
    k(1:1,3) = b+xpos(1:1,3)*(xpos(1:1,1)-c);
end