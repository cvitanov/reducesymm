clear all
hold off
figure(1)
clc

%initial position:
%load ui.mat
%xi = ui;

%xi = [-7.0510; 0.93089; -13.597; -0.50085];
%xi=rand(4,1);
xi = [0.8; 1.2; 0.9; 1.1];
xi= [1.5652;
	 2.5564;
	-4.3934;
	 1.7731];
%xi= [0.0019262; 1.4230754; -0.0010525; 0.0054132];
%xi=[12.1951;6.09756;11.6002;0];

deltat = 0.01;
tfinal = 500;

x = integrator(xi, tfinal, deltat);
save('timeev.mat', 'x', 'tfinal', 'deltat'); % Save the time evolution and parameters

plotflow(0,tfinal,1,2,3);
%view(160,15)

xlabel('u')
ylabel('v')
zlabel('w')

