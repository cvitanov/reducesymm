function  void = integrateMhat(xhati, phii, tinitial, tfinal, deltat, xhatp)
%Integrates the system with the velocity function velocity(x) within the
%slice defined by the template point xhatp according to chaosbook eq 
%10.50.51
%This function reads T-matrix from generator.mat, returns and saves the 
%time evolution of xhat and phi in a file.

load generator.mat %load G

n = floor(tfinal / deltat); % number of time steps

xhat(:,1) = xhati; %x at t=0
phi(:,1) = phii;
xreconstructed(:,1) = LieElement(phi(1), xhat(:,1));

tp = T*xhatp;

vphi = @(xhat) dot(velocity(xhat), tp)/dot(T*xhat, tp);
vhat = @(xhat) velocity(xhat) - vphi(xhat)*(T*xhat);

tol=1e10;

for i = 1:n
	
	if abs(vphi(xhat(:,i))) > tol;
		
		fprintf('encountered a chart border, leaving at:')
		
		t = deltat * (i-1);
		
		break;
		
	end
	
	k1 = deltat*vhat(xhat(:,i));
	k2 = deltat*(vhat(xhat(:,i) + k1/2));
	k3 = deltat*(vhat(xhat(:,i) + k2/2) );
	k4 = deltat*(vhat(xhat(:,i) + k3));
	
	xhat(:,i+1) = xhat(:,i) + k1/6 + k2/3 + k3/3 + k4/6;

	k1phi = deltat*vphi(xhat(:,i));
	k2phi = deltat*(vphi(xhat(:,i) + k1/2));
	k3phi = deltat*(vphi(xhat(:,i) + k2/2) );
	k4phi = deltat*(vphi(xhat(:,i) + k3));
	
	phi(i+1) = phi(i) + k1phi/6 + k2phi/3 + k3phi/3 + k4phi/6;
	
	xreconstructed(:,i+1) = LieElement(phi(i+1), xhat(:,i+1));
	
	if phi(:,i+1) < 0
		phi(:,i+1) = phi(:,i+1) + 2*pi;
	end
	if phi(:,i+1) > 2*pi
		phi(:,i+1) = phi(:,i+1) - 2*pi;
	end
	
end

save('dynamicswithintheslice.mat','xhat','phi','xreconstructed','deltat','tinitial','tfinal')
