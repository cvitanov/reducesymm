clear;
clc;

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

%Parameters variation:
mu1 = -2.8023
a1 =  -1
b1 =  0
c1 =  -4.1122
mu2 = 1
a2 =  -2.6636
b2 =  0
c2 =  1.8854
e2 = 2.1769

%reassign parameters:
pars(1) = mu1;
pars(2) = a1;
pars(3) = b1;
pars(4) = c1;
pars(5) = mu2;
pars(6) = a2;
pars(7) = b2;
pars(8) = c2;  
pars(9) = e2;

save('pars.mat','pars')
