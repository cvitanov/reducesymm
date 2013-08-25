function wandq = wqeq(u,v)

%Function calculates w and q variables from given u, v and parameters

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

A1 = mu1 + a1*u + b1*v;
A2 = mu2 + a2*u + b2*v;

w=(-2*u*A1)/c1;
q=(-e2*w)/(2*A1 + A2);

wandq = [w;q];
