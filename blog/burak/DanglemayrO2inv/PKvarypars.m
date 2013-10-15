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

%mu1=1.768906995148714   
%mu1=1.769   
%mu1=1.9   
mu1=2   
%a1=0.406357259437816  
%a1=0.406 
%a1=0.41 
a1=0.41 
%b1=-1.660768375284920   
%b1=-1.661
%b1=-1.66
b1=-1.7
c1=1.000000000000000
%mu2=-0.675564867234430   
%mu2=-0.676
%mu2=-0.67
%mu2=-0.7
mu2=-1
%a2=0.083129941576966  
%a2=0.0831  
%a2=0.1  
a2=0.16  
%b2=-0.047035211500671  
%b2=-0.0470
b2=-0.05
%b2=-0.04803258
c2=-1.000000000000000
e2=0

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
