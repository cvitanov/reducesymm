function A = StabilityMatrix(x)

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

u = x(1);
v = x(2);
w = x(3);
q = x(4);

A = [2*mu1+4*a1*u+2*b1*v  2*b1*u c1 0;
	 2*a2*v 2*mu2+2*a2*u+4*b2*v c2 0;
	 4*c2*u+4*c1*v+(2*a1+a2)*w 4*c1*u+(2*b1+b2)*w 2*mu1+mu2+(2*a1+a2)*u+(2*b1+b2)*v -e2;
	 (2*a1+a2)*q (2*b1+b2)*q e2 2*mu1+mu2+(2*a1+a2)*u+(2*b1+b2)*v];
