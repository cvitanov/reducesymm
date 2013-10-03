clear;
clc;

%Invariant subspace of Porter-Knobloch system in the invariant polynomial
%basis is given by (0,v,0,0). This program computes corresponding points in
%full state space. Since u = x_1^2 + x_2^2 = 0, x_1=x_2=0.

j = 0;
imax=10;
for v = 0:0.1:14
	
	j=j+1;
	relation = @(y) y(1)^2+y(1)^2-v;
	
	for i = 1:imax
			
		[y, fval, info] = fsolve(relation, randn(2,1));
	
		xinvsubs(:,(j-1)*imax+i) = [0;0;y];
		
	end
end
save xinvsubs.mat xinvsubs
