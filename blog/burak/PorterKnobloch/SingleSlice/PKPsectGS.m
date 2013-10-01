clear all
clc

Rx = @(theta, x) [1 0 0;
				  0 cos(theta) -sin(theta);
				  0 sin(theta) cos(theta)]*x;
				  
Ry = @(theta, x) [cos(theta)  0 sin(theta);
				  0 		  1 0;
				  -sin(theta) 0 cos(theta)]*x;

Rz = @(theta, x) [cos(theta) -sin(theta) 0;
				  sin(theta) cos(theta)  0;
				  0 		 0 			 1]*x;

%Poincare section template:				 

%1st equilibrium
xhatp1=[-4.39840242929180e-01;
	   7.28510213391687e-01;
	   -1.32288747382379e-01];
%Unstable direction of the 1st eq:
vu1 = [0.1755064524073753;
	   -0.9800140886924088;
	   0.0936475900798591];
%A vector perpendicular to xhatp1 & xhatp1+vu1:
nhat1=cross(xhatp1+vu1, xhatp1);
nhat1=nhat1/norm(nhat1); %Normalize	   
	   
ps1 = psectgs(xhatp1, nhat1, 'psect1.mat');
V1 = [vu1 cross(vu1,nhat1) nhat1]; %Similarity matrix.

figure(1);
plotPsectGS('psect1.mat', V1)

xhatp0 = [0;0;0];

vu0 = [0.000000000000000;
	   0.841470984807897;
	   0.540302305868140];

nhat0 = Rx(pi/2, vu0);

ps0 = psectgs(xhatp0, nhat0, 'psect0.mat');
V0 = [cross(vu0,nhat0) vu0 nhat0];

figure(2);
plotPsectGS('psect1.mat', V0)
