clear all
hold off
figure(1)
clc

%EDIT THIS FILE WITH EXTREME CARE
%initial position:
%load xi.mat
%xi = xi;

plot = 1;
xi = [0;
	  0;
	  0;
	  0;
	  1;
	  0;
	  0;
	  0;
	  0;
	  1;
	  0;
	  0;
	  0;
	  0;
	  1;
	  0;
	  0;
	  0;
	  0;
	  1] ; % Starting point close to the relative equilibrium.

%xi = [-0.311014008385794;
	  %-0.311014028442128;
	  %-0.132288747382379;
	  %-0.728510213391687;
	  %1;
	  %0;
	  %0;
	  %0;
	  %0;
	  %1;
	  %0;
	  %0;
	  %0;
	  %0;
	  %1;
	  %0;
	  %0;
	  %0;
	  %0;
	  %1] ; % Starting point close to the relative equilibrium.


deltat = 0.01;
tfinal = 1;
%tfinal = 11.55;

x = integratorvar(xi, tfinal, deltat);

if plot
	plotflowvar(0, tfinal, 1,2,3);
	xlabel('x_1')
	ylabel('x_2')
	zlabel('y_1')
	view(120,30)
end
