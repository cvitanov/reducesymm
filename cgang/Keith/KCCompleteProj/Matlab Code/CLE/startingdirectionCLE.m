%Written by K. Carroll 5/1/2012
function [ v1 ] = startingdirectionCLE(xs)
    b = 8/3;
    p1 = 28;
    e = 1/10;
    sig = 10;
    d = 1+e^2/((sig+1)^2);
    r = sqrt(b*(p1-d));
    r2 = sqrt(d)*r;
    z = p1-d
    theta = -acos(1/sqrt(d));
    
    A = [-sig, 0, sig, 0, 0,
        0, -sig,  0, sig, 0,
        p1-xs(1:1,5), 0, -1, -e, -xs(1:1,1),
        0, p1-xs(1:1,5), e, -1, -xs(1:1,2),
        xs(1:1,3), xs(1:1,4), xs(1:1,1), xs(1:1,2), -b];
    [V,D ] =eig(A);
    for i=1:5
        if real(D(i:i,i:i)>0)
            v1 = real(V(1:5,i));
            break;
        end
    end
    
end