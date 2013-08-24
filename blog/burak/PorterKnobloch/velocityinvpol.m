function v = velocityinvpol(x)
%Porter - Knobloch velocity function

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

v = [2*mu1*x(1) + 2*a1*x(1)^2 + 2*b1*x(1)*x(2) + c1*x(3);
	 2*mu2*x(2) + 2*a2*x(1)*x(2) + 2*b2*x(2)^2 + c2*x(3);
	(2*mu1+mu2)*x(3) + (2*a1+a2)*x(1)*x(3) + (2*b1+b2)*x(2)*x(3) + 4*c1*x(1)*x(2)+2*c2*x(1)^2-e2*x(4);
	(2*mu1 + mu2)*x(4) + (2*a1 + a2)*x(1)*x(4) + (2*b1 + b2)*x(2)*x(4) + e2*x(3)];
