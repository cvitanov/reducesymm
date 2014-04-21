function [th_p, ph_p, grot] = hexfunc(th, ph)
% this function uses the un-naturla variables

d = sqrt(3) * sin(th + ph) / (2 * cos(th));
R = sqrt(3)/2 * 1/1.15;

grot = 0;
if abs(d) > R
    % normal hex mapping here
    [tanth_p, tanph_p, grot] = freehexfunc(th, ph);
    th_p = atan(tanth_p);
    ph_p = atan(tanph_p);
    return;
end


varph = asin(d/R) - (th+ph);
angchg = scatteringfunc(d, R);
ph_1 = ph + angchg;
if d >= 0
    th_1 = 2*(th+ph+varph) - (th+varph) + angchg;
else
    th_1 = 2*(th+ph+varph) - (th+varph) + angchg + 2*pi;
end

[th_1, ph_1, grot1] = hexgrpaction(th_1, ph_1);
[th_p, ph_p, grot2] = circle2hex(th_1, ph_1, R);
grot = grot1 + grot2;


end

