clear all
hold off
figure(1)
clc

%initial position:
%load xi.mat
%xi = xi;

plot = 1;
xi = [-0.310563473904;
	  -0.310563473904;
	  -0.0419063666849;
	  -0.38981313744] ; % Starting point close to the relative equilibrium.

deltat = -0.01;
tfinal = -0.5;

x = integrator(xi, tfinal, deltat);

save('timeev.mat', 'x', 'tfinal', 'deltat'); % Save the time evolution and parameters

if plot
	plotflow(0, tfinal, 1,2,3);
	xlabel('x_1')
	ylabel('x_2')
	zlabel('y_1')
	view(120,30)
end
