%Written by K. Carroll 5/1/2012
function [ V ] = crossover(x1, x2, n, pos)
    d1 = x1-pos;
    d2 = x2-pos;
    dot1 = dot(d1,n);
    dot2 = dot(d2,n);
    if (dot1<0) && (dot2>0)
        V = 1;
    else
        V = 0;
    end
end