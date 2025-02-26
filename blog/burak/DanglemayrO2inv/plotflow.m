function void = plotflow(ti,tf,i,j,k);

hold off;
clc
figure(1); % Plot on figure(1)
load timeev.mat; % Load x(t), deltat

ni = floor(ti/deltat)+1;
nf = floor(tf/deltat)+1;

if i && j && k
	plot3(x(i, ni:nf), x(j, ni:nf), x(k, ni:nf));
else
	plot3(x(1, ni:nf), x(2, ni:nf), x(3, ni:nf));
end
view(120,30);
box off;
void = 1;
