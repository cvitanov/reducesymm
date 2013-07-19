function dx = EOMPolar(t,x)
dx = zeros(3,1);

mu1 = -0.38;
mu2 = 0.38;
c1 = -1;
c2 = 1;
a1 = -1.31;
a2 = -2.6;
b1 = 1.504;
b2 = 0.22;
e2 = 1.61;

r1 = x(1);
r2 = x(2);
phi = x(3);

dx(1) = r1*(a1*r1^2 + b1*r2^2 + c1*cos(phi)*r2 + mu1);
dx(2) = a2*r1^2*r2 + c2*cos(phi)*r1^2 + b2*r2^3 + mu2*r2;
dx(3) = e2 - 2*c1*r2*sin(phi) - (c2*r1^2*sin(phi))/r2;