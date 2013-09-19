function void = plotreducedflow(ti,tf,i,j,k,options);

load timeev.mat;
clear x;
deltat = deltat;
load movingframes.mat; % Load x(t), deltat

ni = floor(ti/deltat)+1;
nf = floor(tf/deltat)+1;

plot3(xhat(i, ni:nf), xhat(j, ni:nf), xhat(k, ni:nf));

view(120,30);

void = 1;
