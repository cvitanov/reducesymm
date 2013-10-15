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


deltat = 0.01;
tfinal = 100;

x = integrator(xi, tfinal, deltat);

plotflow(0,tfinal,1,2,3);
%view(160,15)

xlabel('u')
ylabel('v')
zlabel('w')
