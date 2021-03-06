function void = MovingFrames2(xhatp)
%Reduces the symmetry 'LieElement' of the system using post processing method
%Saves resulting phi and xhat in movingframes.mat output file

load('timeev.mat'); % Load x, tfinal, deltat integrated by RK4
load('generator.mat'); %Load the Lie element generator

%Template tangent
tp = T * xhatp;
%Slice condition:
SliceCondition = @(phi, x) dot(tp, LieElement(-phi, x)); %eq 10.44

for i = 1:size(x,2)

i;

deltaphi = 0.1; %starting phi-step

	j=0;

	for phi=0:deltaphi:2*pi
		
		if SliceCondition(phi, x(:,i)) < 0 && SliceCondition(phi+deltaphi, x(:,i))>0
			
			j++;
			phidummy(j) = phi;
		    xcandidate(:,j) = LieElement(-phi,x(:,i));
		    xdistsquare(j) = sum((xhatp-xcandidate(:,j)).^2);
		    
		end
		
	end
	
	clear xcandidate;
	
	jstar = find(xdistsquare==min(xdistsquare));
	
	clear xdistsquare;
	
	phi = phidummy(jstar);
	
	clear phidummy;
	
	SliceCondFixedx = @(phi) dot(tp, LieElement(-phi, x(:,i)));
	
	[phii, fval, info] = fsolve(SliceCondFixedx, phi);

	phitau(i) = phii;
	xhat(:,i) = LieElement(-phii,x(:,i));

end

save('movingframes.mat','phitau','xhat')
