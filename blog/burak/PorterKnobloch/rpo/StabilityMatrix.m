function stabmax = StabilityMatrix(x)

%load parameters:
load pars.mat
%assign parameters;
mu1 = pars(1); 
a1 =  pars(2);
b1 =  pars(3);
c1 =  pars(4);
mu2 = pars(5);
a2 =  pars(6);
b2 =  pars(7);
c2 =  pars(8);
e2 =  pars(9);

x1 = x(1);
x2 = x(2);
y1 = x(3);
y2 = x(4);

stabmax = [3*x1^2*a1 + x2^2*a1 + y1^2*b1 + y2^2*b1 + y1*c1 + mu1  2*x1*x2*a1 + y2*c1  2*x1*y1*b1 + x1*c1  2*x1*y2*b1 + x2*c1;
		   2*x1*x2*a1 + y2*c1  x1^2*a1 + 3*x2^2*a1 + y1^2*b1 + y2^2*b1 - y1*c1 + mu1  2*x2*y1*b1 - x2*c1  2*x2*y2*b1 + x1*c1;
		   2*x1*y1*a2 + 2*x1*c2  2*x2*y1*a2 - 2*x2*c2  x1^2*a2 + x2^2*a2 + 3*y1^2*b2 + y2^2*b2 + mu2  2*y1*y2*b2+e2;
		   2*x1*y2*a2 + 2*x2*c2  2*x2*y2*a2 + 2*x1*c2  2*y1*y2*b2 - e2   x1^2*a2 + x2^2*a2 + y1^2*b2 + 3*y2^2*b2  + mu2];
