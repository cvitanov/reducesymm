clear
clc;

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

eq = [0;-mu2/b2;0;0];

lambda = eig(StabilityMatrix(eq))
