clear
clc
hold off

%Load the geometry package for plotting vectors:
pkg load geometry;

%Template point:
xhatp=[1;1;0;0];

load gramschmidt.mat; %Load the flow on the Gram-Schmidt Coordinates

%Place a red dot on the unstable relative equilibrium
plot3(xhatGS(1,1), xhatGS(2,1), xhatGS(3,1), 'r.', 'markersize', 8)

hold on

%Plot the spiral-out from the relative equilibrium green:
plotgramschmidt(0.01,250,1,2,3,'g');

%Plot the flow on the Gram-Schmidt coordinates:
plotgramschmidt(250,500,1,2,3,'b');

%Place a black dot on the unstable equilibrium at the origin
plot3(0,0,0,'k.', 'markersize', 8)


%Load Jacobian, Jacobian eigenvectors and eigenvalues calculated for
%the unstable relative equilibrium at the initial point of the flow:

load twandjac.mat;

%Calculate the unstable direction:

vu = L(3,3)*V(:,3)+L(4,4)*V(:,4);

%Remove the component on the direction of the group tangent:
load generator.mat;
t=T*xhatp;
t=t/norm(t);
vu = vu - (t'*vu)*t;
%Project it onto the Gram-Schmidt Coordinates:
vuGS = findGS(vu);
vuGS=vuGS/norm(vuGS);
hvu=drawVector3d(xhatGS(:,1)', vuGS');
set(hvu, 'linewidth', 2)
set(hvu,'color','r')

%Calculate the unstable direction of the origin:
load origin.mat;
vu0 = exp(L0(1,1))*V0(:,1) + exp(L0(2,2))*V0(:,2);
vu0GS=findGS(vu0);
vu0GS=vu0GS/norm(vu0GS);
hvu0=drawVector3d([0 0 0], vu0GS');
set(hvu0, 'linewidth', 2)
set(hvu0,'color','k')
view(280,30);

ah = gca();%Generate a handler for the current Axis object
set(ah,'box','off', 		%Remove the box.
	   'xtick',[-2.5:0.5:0.5],
	   'xticklabel',{'','-2','','-1','','0',''},
	   'ytick',[-1:5],
	   'yticklabel',{'','0','','2','','4',''},
	   'ztick',[-1.5:0.5:1.5],
	   'zticklabel',{'','-1','','0','','1',''}) 

print('unstablevectors.eps','-S650,450',
'-depsc2',
'-F:Helvetica:8'
)
