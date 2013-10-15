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

view(25,25);
ah = gca();%Generate a handler for the current Axis object
set(ah,'box','off') 		%Remove the box.
void = 1;
