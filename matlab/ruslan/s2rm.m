function s1 = s2rm(s, npts )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

s1=zeros(npts(end),2);
for i=1:npts(end),
    s1(i,:) = [s(i), s(npts(1)+i)];
end

end

