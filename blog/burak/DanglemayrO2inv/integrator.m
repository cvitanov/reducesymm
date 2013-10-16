function x = integrator(xi, tfinal, deltat)
%Runge-Kutta 4 (ChaosBook exercise 2.6) for the system with velocity 
%function velocity(x). 
%xi = initial point for integration
%t = duration of integration
%deltat = length of timesteps

n = floor(tfinal / deltat); % number of time steps

x(:,1) = xi; %x at t=0

for i = 1:n

	k1 = deltat*velocity(x(:,i));
	k2 = deltat*(velocity(x(:,i) + k1/2));
	k3 = deltat*(velocity(x(:,i) + k2/2) );
	k4 = deltat*(velocity(x(:,i) + k3));
	
	x(:,i+1) = x(:,i) + k1/6 + k2/3 + k3/3 + k4/6;

end
