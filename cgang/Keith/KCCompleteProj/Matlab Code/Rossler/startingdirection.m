%Written by K. Carroll 5/1/2012
function [ v1 ] = startingdirection(xs)
    a = 0.2;
    b = 0.2;
    c = 5.7;
    A = [0, -1, -1,
        1, a,  0,
        xs(1:1,3), 0, (xs(1:1,1)-c)];
    [V, D ] =eig(A)
    for i=1:3
        if real(D(i:i,i:i)>0)
            v1 = real(V(1:3,i));
            break;
        end
    end
    
end