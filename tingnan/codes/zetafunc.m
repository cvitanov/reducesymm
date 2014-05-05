function [ zetaval ] = zetafunc(np, tp)

zetaval = 1;
for ii = 1:np 
    tmpmat = combntns(tp, ii);
    tmpvec = prod(tmpmat, 2);
    sumval = sum(tmpvec)*(-1)^ii;
    zetaval = zetaval + sumval;
end
end

