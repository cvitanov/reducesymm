clear
clc

load pars.mat
mu1 = pars(1);
q1 = pars(2);
mu2 = pars(3);
q2 = pars(4);
e2 = pars(5);

equilibria(:,1) = [0;
				   0;
				   0;
				   0];
				 
equilibria(:,2) = [-((2*mu1*mu2*(e2^2 + (2*mu1 + mu2)^2))/((4*mu1^2 - mu2^2)*q1*q2));
				   -((mu1^2*(e2^2+(2*mu1 + mu2)^2))/((4*mu1^2 - mu2^2)*q1^2));
				   (4*e2*mu1^2*mu2*(e2^2+(2*mu1+mu2)^2))/((2*mu1-mu2)*(2*mu1+mu2)^2*q1^2*q2);
				   (4*mu1^2*mu2*(e2^2+(2*mu1+mu2)^2))/((4*mu1^2+mu2^2)*q1^2*q2)]

for i = 1 : size(equilibria,2)				   
	A = StabilityMatrix(equilibria(:,i));
	lambda(:,i) = eig(A);
end

lambda = lambda


