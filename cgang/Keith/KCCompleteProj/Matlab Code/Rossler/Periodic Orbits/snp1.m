function [sn1] = snp1(sn, c)
    s1 = size(c);
    sn1 = 0;
    for i=1:s1(1:1,1)
        sn1 = sn1 + c(i:i,1)*(sn)^(i-1);
    end
end