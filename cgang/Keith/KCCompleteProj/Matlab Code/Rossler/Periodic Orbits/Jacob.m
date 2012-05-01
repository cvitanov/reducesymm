% Written by K. Carroll 5/1/2012
function [ M ] = Jacob(M1, x1, deltat)
    a = 0.2;
    b = 0.2;
    c = 5.7;

    A = [0, -1, -1,
        1, a, 0
        x1(1:1,3), 0, x1(1:1,1)-c];
    M2 = eye(3) + A*deltat;
    M = M2*M1;
end