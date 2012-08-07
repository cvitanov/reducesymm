function dx = EOMCartesian(t,x)
dx = zeros(4,1);

mu1 = -0.38;
mu2 = 0.38;
c1 = -1;
c2 = 1;
a1 = -1.31;
a2 = -2.6;
b1 = 1.504;
b2 = 0.22;
e2 = 1.61;

x1 = x(1);
x2 = x(2);
y1 = x(3);
y2 = x(4);

dx(1) = a1*x1^3 + a1*x1*x2^2 + b1*x1*y1^2 + c1*x1*y1 + b1*x1*y2^2 + mu1*x1 + c1*x2*y2;
dx(2) = a1*x1^2*x2 + c1*x1*y2 + a1*x2^3 + b1*x2*y1^2 - c1*x2*y1 + b1*x2*y2^2 + mu1*x2;
dx(3) = a2*x1^2*y1 + c2*x1^2 + a2*x2^2*y1 - c2*x2^2 + b2*y1^3 + b2*y1*y2^2 + mu2*y1 + e2*y2;
dx(4) = a2*x1^2*y2 + 2*c2*x1*x2 + a2*x2^2*y2 + b2*y1^2*y2 - e2*y1 + b2*y2^3 + mu2*y2;