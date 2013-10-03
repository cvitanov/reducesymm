function void = plotgramschmidt(ti,tf,i,j,k,options);

load timeev.mat;
clear x;
deltat = deltat;
load gramschmidt.mat; % Load x(t), deltat

ni = floor(ti/deltat)+1;
nf = floor(tf/deltat)+1;

plot3(xhatGS(i,ni:nf), xhatGS(j,ni:nf), xhatGS(k,ni:nf), options)
xlabel('x_{GS}')
ylabel('y_{GS,1}')
zlabel('y_{GS,2}')
title('Flow on Gram-Schmidt Coordinates')

view(285,25);

void = 1;
