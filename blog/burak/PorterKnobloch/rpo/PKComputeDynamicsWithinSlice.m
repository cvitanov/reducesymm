clear
hold off
clc

xhatp=[1;1;0;0]; %Slice template
movingframes=0;

if movingframes

	load movingframes.mat;
	xhati=xhat(:,1); %Initial point
	phii=phitau(1); %Initial angle
	clear xhat;
	%clear phitau;
	
else
	
	xhati = [-0.583132295203164;
			 -0.583132295203164;
			 -0.121331636628869;
			 -2.761084984756857];
	  
	phii=0; %Initial angle	
	
end

tinitial = 0;
tfinal = 150;
deltat = 0.01;

integrateMhat(xhati, phii, tinitial, tfinal, deltat, xhatp);
