clear
hold off
clc

xhatp=[1;1;0;0]; %Slice template
movingframes=1;

if movingframes

	load movingframes.mat;
	xhati=xhat(:,1); %Initial point
	phii=phitau(1); %Initial angle
	clear xhat;
	%clear phitau;
	
else
	
	xhati=xhatp; %Initial point
	phii=0; %Initial angle	
	
end

tinitial = 0;
tfinal = 100;
deltat = 0.01;

integrateMhat(xhati, phii, tinitial, tfinal, deltat, xhatp);
