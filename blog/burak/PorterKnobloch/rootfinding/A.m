function stabmax = A(x)

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

stabmax = [3*x(1)^2*a1 + x(2)^2*a1 + x(3)^2*b1 + x(4)^2*b1 + x(3)*c1 + mu1, 2*x(1)*x(2)*a1 + x(4)*c1, 2*x(1)*x(3)*b1 + x(1)*c1, 2*x(1)*x(4)*b1 + x(2)*c1;
		   2*x(1)*x(2)*a1 + x(4)*c1, x(1)^2*a1 + 3*x(2)^2*a1 + x(3)^2*b1 + x(4)^2*b1 - x(3)*c1 + mu1, 2*x(2)*x(3)*b1 - x(2)*c1, 2*x(2)*x(4)*b1 + x(1)*c1;
		   2*x(1)*x(3)*a2 + 2*x(1)*c2, 2*x(2)*x(3)*a2 - 2*x(2)*c2, x(1)^2*a2 + x(2)^2*a2 + 3*x(3)^2*b2 + x(4)^2*b2 + mu2, 2*x(3)*x(4)*b2+e2;
		   2*x(1)*x(4)*a2 + 2*x(2)*c2, 2*x(2)*x(4)*a2 + 2*x(1)*c2, 2*x(3)*x(4)*b2 - e2 , x(1)^2*a2 + x(2)^2*a2 + x(3)^2*b2 + 3*x(4)^2*b2  + mu2];
