%Written by K. Carroll 5/1/2012
function [ x3 ] = crossoverpos(x1, x2, n, pos)
    diff = x2-x1;
    planedot = dot(n,pos);
    pos1dot = dot(n,x1);
    diffdot = dot(n,diff);
    t = (planedot-pos1dot)/diffdot;
    x3 = x1+diff*t;
end