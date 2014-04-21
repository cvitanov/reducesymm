function [ x ] = hexmapping(x0)
% We will use the natural variables, the position ll on the hexboundary, 
% and the tangent component of velocity, sin(ph), to continue our
% calculation. The function are expecting to take a 2-dimensional vector.

th = x0(1);
ph = x0(2);
[th, ph, gaction] = hexfunc(th, ph);
% if th <= 0
%     th = -th;
%     ph = -ph;
% end
% 
% if th == pi/6 && ph>0
%    ph = -ph;
% end

x = zeros(2,1);
x(1) = th;
x(2) = ph;

end

