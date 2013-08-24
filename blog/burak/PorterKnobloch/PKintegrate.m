clear all
hold off
figure(1)
clc

%initial position:
load xi.mat
xi = xi;

%xi = [9.36; 5.64; -44.46; 0.10] ;

deltat = 0.01;
tfinal = 100;

x = integrator(xi, tfinal, deltat);

plot3(x(1,:), x(3,:), x(4,:));
xlabel('x_1')
ylabel('y_1')
zlabel('y_2')
view(120,30)
