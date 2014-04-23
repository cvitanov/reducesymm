function [ x1 ] = hexmappingelem( x0 )
%HEXMAPPINGELEM Summary of this function goes here
%   Detailed explanation goes here
ith = x0(1);
iph = x0(2);
R = sqrt(3)/2*1/1.15;
[th1, ph1, grot1] = hexgrpaction(ith, iph+ith);
[thh, phh, grot2] = circle2hex(th1, ph1, R);
[thh, phh, groth] = hexgrpaction(thh, phh);
grot = grot1 + grot2 + groth;
d = sqrt(3) * sin(thh + phh) / (2 * cos(thh));
while abs(d) > R
    % normal hex mapping here until the things hit the circle
    [tanth, tanph, tmpgrot] = freehexfunc(thh, phh);
    thh = atan(tanth);
    phh = atan(tanph);
    d = sqrt(3) * sin(thh + phh) / (2 * cos(thh));
    grot = grot + tmpgrot;
end
varph = asin(d/R) - (thh+phh);
angchg = scatteringfunc(d, R);
ph2 = phh + angchg;
if d >= 0
    th2 = 2*(thh+phh+varph) - (thh+varph) + angchg;
else
    th2 = 2*(thh+phh+varph) - (thh+varph) + angchg + 2*pi;
end


fth = mod(th2+grot*pi/3, 2*pi);
fph = mod(ph2-th2, 2*pi);

if fph >= pi/2*3
    fph = fph - 2*pi;
end


x1 = [fth;fph];
% if  fph > pi/2 || abs(fth) > pi/6
%     ith
%     iph
%     error('bail out')
% end


end

