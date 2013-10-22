clear all
hold off
figure(1)
clc

%initial position:
%load xi.mat
%xi = xi;

plot = 1;
%xi = [0.181612;
	 %-0.512500;
	  %1.939339;
	 %-0.010343] ; % Starting point close to the relative equilibrium.
xi = [-3.77393;
	 -1.55696;
	  -2.35925;
	 1.66352] ; % Starting point close to the relative equilibrium.

deltat = 0.01;
tfinal = 100;

x = integrator(xi, tfinal, deltat);

save('timeev.mat', 'x', 'tfinal', 'deltat'); % Save the time evolution and parameters

if plot
	plotflow(0, tfinal, 1,2,3);
	xlabel('x_1')
	ylabel('x_2')
	zlabel('y_1')
	view(120,30)
end
