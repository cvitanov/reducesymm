function ps=psectgs(xhatp, nhat, filename, direction)

load gramschmidt.mat;
load timeev.mat;

%Plane condition:

%direction=-1; %Now a function input

U = @(x) dot((x-xhatp),nhat)*direction;

n = size(xhatGS,2);
tol= 1e-6; %Tolerance for 'zero'

ps = [xhatp;0]; %Assign a dummy value to ps vector to be able to fill afterwards

for i = 1:n-1
	%If plane equation change sign from negative to positive, look for the intersection:
	if U(xhatGS(:,i)) < 0 && U(xhatGS(:,i+1)) > 0 
		%Integration until the plane condition is satisfied:
			
		xdummyGS = xhatGS(:,i);
		xdummyfull = x(:,i);

		DotProduct = U(xdummyGS); 
		deltatadaptive = deltat/2;
		pstime = (i-1)*deltat;
		j=1	;
		while abs(DotProduct) > tol
			j=j+1;
			%Adaptive integration: 
			
			xint=integrator(xdummyfull, deltatadaptive, deltatadaptive);
			
			
			pstime = pstime + deltatadaptive;
			
			xdummyfull = xint(:,2);
			
			xdummyhat = findxhat(xdummyfull,[1;1;0;0]);
			xdummyhat = xdummyhat(1:size(xdummyhat)-1);
			
			xdummyGS = findGS(xdummyhat);

			DotProduct = U(xdummyGS);
			
			% If the Poincare section is passed get back in integration
			
			if DotProduct > 0 % If the Poincare section is passed
				
				xdummyfull = xint(:,1);
				pstime = pstime - deltatadaptive;
				deltatadaptive = deltatadaptive/2;
				
				%xdummyhat = findxhat(xdummyfull,[1;1;0;0])
				%xdummyhat = xdummyhat(1:size(xdummyhat)-1)
				%xdummy = findGS(xdummyhat)
				%DotProduct = U(xdummy)
				
			end
		
			if j>100
				fprintf('not converging, exiting...\n')
				break;
			end
		end
		
		ps(:,size(ps,2)+1) = [xdummyGS;
							  pstime];
		%ps(:,size(ps,2)+1) = [xdummyGS;
							  %pstime];
	end
	
end

ps = ps(:, 2:size(ps,2));
save(filename, 'ps','xhatp')
