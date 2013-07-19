% Written by K. Carroll 5/1/2012
function [A] = testifsame(v1,v2)
    s = size(v1);
    colnum = s(1:1,2);
    M = zeros(colnum, colnum);
    for i=1:colnum
        for j=1:colnum
            addv = j+i-1;
            addvmod = mod(addv, colnum)+1;
            M(i:i,j:j) = v1(1, addvmod:addvmod); 
        end
    end
    A = 0;
    for i=1:colnum
        if (M(i,:)==v2)
            A =1;
        end
    end
end