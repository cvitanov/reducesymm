clear all
hold off
figure(1)
clc

%initial position:
%load xi.mat
%xi = xi;

plot = 1;
xi = [0.12712; -0.42107; 0.68001; 0.29294] ; % Starting point close to the relative equilibrium.

deltat = 0.01;
tfinal = 1000;

x = integrator(xi, tfinal, deltat);

if plot
	plotflow(0, tfinal, 1,2,3);
	xlabel('x_1')
	ylabel('x_2')
	zlabel('y_1')
	view(120,30)
end
