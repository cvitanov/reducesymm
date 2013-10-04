clear all
hold off
figure(1)
clc

%initial position:
%load xi.mat
%xi = xi;

plot = 1;
  
xi = [-0.583132450887042;
	  -0.583132450887042;
	  -0.121331644065296;
	  -2.761085488025118] ; % Starting point close to the relative periodic orbit.

deltat = 0.01;
tfinal = 150;

x = integrator(xi, tfinal, deltat);

save('timeev.mat', 'x', 'tfinal', 'deltat'); % Save the time evolution and parameters

if plot
	plotflow(0, tfinal, 1,2,3);
	xlabel('$x_1$')
	ylabel('$x_2$')
	zlabel('$y_1$')
	view(160,20)
	
	print -depslatexstandalone flowfullssp.tex
	
end
