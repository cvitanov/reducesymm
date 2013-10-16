function ps = poincare(xhatp, nhat)
%Poincare section calculator using method described in ChaosBook sect 3.2 
%xname=string to tell the name of file keeping the time evolved x values
%xhatp = a point on the poincare section
%n = normal of the poincare section

load('timeev.mat'); % Load x, tfinal, deltat integrated by RK4

%Plane equation for the Poincare section (ChaosBook eq 3.6):

U = @(x) dot((x-xhatp),nhat);

n=size(x,2);
N=10; % Number of divisions to look for intersection
tol=eps; %Tolerance for 'zero'
ps = [xhatp;0]; %Assign a dummy value to ps vector to be able to fill afterwards

for i = 1:n-1
	%If plane equation change sign from negative to positive, look for the intersection:
	if U(x(:,i)) > 0 && U(x(:,i+1)) < 0 
	
		%Choose the coordinate with the maximum speed at the instance
		%as the new time coordinate for the integration to avoid 
		%singular velocities
		
		vx = velocity(x(:,i)); % velocity at this point
		absvx = abs(vx); % absolut value of the velocity
		j = find(absvx == max(absvx)); 
			
		%Integration until the plane condition is satisfied:
			
		xdummy = x(:,i);
		deltaxj = (x(j,i+1) - x(j,i))/N;
		xj = xdummy(j);
			
		%Transformed state space vector:
		xnew = xdummy;
		xnew(j) = (i-1)*deltat; %jth variable is time!
		k = 0;
		l=0;		
		pstime = xnew(j);
		
		DotProduct = U(xdummy); 
		
		while DotProduct > tol
		
			%Adaptive Euler integration: 
			xnew = xnew + deltaxj * velocityPS(xdummy,j);
			xj = xj + deltaxj; % update xj
			
			xdummy2 = xdummy; % Second dummy variable for backwards integration
			xdummy = xnew;
			xdummy(j) = xj;
			pstime = xnew(j);
			
			DotProduct = U(xdummy);
			
			% If the Poincare section is passed get back in integration
			
			if DotProduct < 0 % If the Poincare section is passed
			
				xnew = xnew - deltaxj * velocityPS(xdummy2,j);
				xj = xj - deltaxj; % update xj
				deltaxj = deltaxj / 2; 
			
			end
		
		end
		
		ps(:,size(ps,2)+1) = [xdummy;
							  pstime];
		
	end
	
end

ps = ps(:, 2:size(ps,2));
