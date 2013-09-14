clear all;
clc;

%pargen();
random = 0;
equilibria = [0;0;0;0]; % First equilibrium is the origin. 
equilibriax = [0;0;0;0]; % First equilibrium is the origin. 
lambda = eig(Ainvpol(0,0,0,0));
lambdax = eig(A([0;0;0;0]));

for ui=0:0.1:1
	for vi=0:0.1:1
	
		if random
		
			uui = 10^(1+rand) * rand;
			vvi = 10^(1+rand) * rand;
			
		else
		
			uui = ui;
			vvi = vi;
			
		end
	
		[uv, fval, info] = fsolve(@fg, [uui;vvi]); % find an equilibrium point
		
		if info == 1 && uv(1)>0 && uv(2)>0
			
			newroot = 1;
			
			for k = 1:size(equilibria,2) % check if the found root is the same with the previous ones
			
				if norm(uv-equilibria(1:2,k)) < 1e-1
			
					newroot = 0;
					
				end
				
			end
			
			if newroot %A new root! Calculate eigenvalues:
					
				u = uv(1);
				v = uv(2);
				w = w = wqeq(u,v)(1);
				q = wqeq(u,v)(2);
				
				uroot=[u;v;w;q];
				save uroot.mat uroot
				
				equilibria(:, size(equilibria,2)+1) = uroot;
				lambda(:,size(lambda,2)+1) = eig(Ainvpol(u,v,w,q));
				
				[x, fval, info] = fsolve(@fxy, randn(4,1));
				equilibriax(:, size(equilibriax,2)+1) = x;
				lambdax(:,size(lambdax,2)+1) = eig(A(x));
				
			end					
		end	
	end
end

%syzeq=syzygy(u,v,w,q)

%ui=[u;v;w;q];
%save ui.mat ui

%[x, fval, info] = fsolve(@fxy, randn(4,1))
%save xi.mat x

%lambdax = eig(A(x))
