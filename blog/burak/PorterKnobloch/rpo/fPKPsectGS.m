function void = fPKPsectGS()

Rx = @(theta, x) [1 0 0;
				  0 cos(theta) -sin(theta);
				  0 sin(theta) cos(theta)]*x;
				  
Ry = @(theta, x) [cos(theta)  0 sin(theta);
				  0 		  1 0;
				  -sin(theta) 0 cos(theta)]*x;

Rz = @(theta, x) [cos(theta) -sin(theta) 0;
				  sin(theta) cos(theta)  0;
				  0 		 0 			 1]*x;

%Poincare section template:				 

%1st equilibrium
xhatp1=[-4.39840242929180e-01;
	   7.28510213391687e-01;
	   -1.32288747382379e-01];
%Unstable direction of the 1st eq:
vu1 = [0.1755064524073753;
	   -0.9800140886924088;
	   0.0936475900798591];
%A vector perpendicular to xhatp1 & xhatp1+vu1:
%nhat1=cross(xhatp1+vu1, xhatp1);

nhat1 = [8;-1;1];
%nhat1 = [1;-1;1]
nhat1 = nhat1 - (dot(vu1,nhat1)/dot(vu1,vu1)) * vu1;

nhat1=nhat1/norm(nhat1); %Normalize	   
	   
ps1 = psectgs(xhatp1, nhat1, 'psect1.mat', -1);
V1 = [-vu1 cross(vu1,nhat1) nhat1]; %Similarity matrix.

figure(1);
plotPsectGS('psect1.mat', V1)
xlabel('$- \hat{v}_{u,TW,GS}$')
ylabel('$\hat{v}_{u,TW,GS} \times \hat{n}_{TW,GS}$')

print -depslatexstandalone PSECT1.tex

figure(3);
hold off;
returnmap('psect1.mat', V1)

print -depslatexstandalone RETMAP1.tex
