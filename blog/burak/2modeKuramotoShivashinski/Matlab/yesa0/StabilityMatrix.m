function stabmax = StabilityMatrix(x)

%load parameters:
load pars.mat
%assign parameters;
mu1 = pars(1);
q1 = pars(2);
mu2 = pars(3);
q2 = pars(4);
e2 = pars(5);
a0 = pars(6);

u = x(1);
v = x(2);
w = x(3);
q = x(4);

stabmax = [2*mu1 		 0  	 0				q1; 
		   0 	 		 2*mu2	 0				q2/2; 
		   0 	 		 0       2*mu1+mu2 		-e2+2*a0*q1-a0*q2;
		   2*q2*u-4*q1*v -4*q1*u e2-a0*q1+a0*q2 2*mu1+mu2];
