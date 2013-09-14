function Jac = Jacobian(ti, tf)
%Computes Jacobian according to ChaosBook eq 4.37

%Load the time evolution (x), total simulation time (tfinal) and 
%length of the timesteps (deltat):

load('timeev.mat')

if tf > tfinal
	fprintf('Error: tf > tfinal. Simulate again until tf.\n');
	return;
end

ni = floor(ti / deltat)+1 %starting step
nf = floor(tf / deltat)+1 %final step

Jac = expm(deltat*StabilityMatrix(x(:,ni)));
%Jac = zeros(size(StabilityMatrix(x(:,ni)),1));

for i = ni+1:nf-1
%for i = ni:nf
	
	i
	Jac =  expm(deltat*StabilityMatrix(x(:,i))) * Jac;
	
%	Jac = Jac + deltat * StabilityMatrix(x(:,i));
	
end

%Jac = expm(Jac);
