clear all
hold off
figure(1)
clc

%initial position:
load ui.mat
xi = ui;

%xi = [9.36; 5.64; -44.46; 0.10] ;

deltat = 0.01;
tfinal = 100;

x = integrator(xi, tfinal, deltat);

plotflow(0,tfinal,1,2,3);

xlabel('u')
ylabel('v')
zlabel('w')
