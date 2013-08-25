function eq18 = fg(uv)

%Equations 18 of 2modes paper

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

u = uv(1);
v = uv(2);

utilda = c2 * u;
vtilda = c1 * v;
a1tilda = a1/c2;
b1tilda = b1/c1;
a2tilda = a2/c2;
b2tilda = b2/c1;

A1 = mu1 + a1tilda * utilda + b1tilda * vtilda;
A2 = mu2 + a2tilda * utilda + b2tilda * vtilda;

eq18 = [utilda*A1 - vtilda*A2;
		(A1^2 - c1*vtilda)*(utilda + 2*vtilda)^2 + e2^2 * vtilda^2];
