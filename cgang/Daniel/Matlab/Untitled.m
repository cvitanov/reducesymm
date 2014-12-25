clear all

T = zeros(4,4);

T(1,2) = 1;
T(2,1) = -1;
T(3,4) = 2;
T(4,3) = -2;

xprime = [1 0 0 0]'

syms mu1 mu2
syms x1 y1 x2 y2 a1 a2 b1 b2 c1 c2 e1 e2
params = [mu1 mu2 a1 a2 b1 b2 c1 c2, e2, 1, 1,1]

tgprime = T*xprime(:)

x = transpose([x1 y1 x2 y2])
tg = T*x

v(1) = a1*x1^3 + b1*x1*x2^2 + c1*x1*x2 + a1*x1*y1^2 + b1*x1*y2^2 + mu1*x1 + c1*y1*y2
v(2) = a1*x1^2*y1 + c1*x1*y2 + b1*x2^2*y1 - c1*x2*y1 + a1*y1^3 + b1*y1*y2^2 + mu1*y1
v(3) = c2*(x1^2 - y1^2) + e2*y2 + mu2*x2 + x2*(a2*x1^2 + a2*y1^2) + 2*b2*x2*y2^2 + b2*x2*(x2^2 - y2^2)
v(4) = a2*x1^2*y2 + 2*c2*x1*y1 + b2*x2^2*y2 - e2*x2 + a2*y1^2*y2 + b2*y2^3 + mu2*y2

phidot =  ((v)*tgprime)/(transpose(tg)*tgprime)

vhat = v - phidot*transpose(tg)