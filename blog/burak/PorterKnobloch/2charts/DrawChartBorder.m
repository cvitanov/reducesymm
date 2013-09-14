clear all;
clc;

load('generator.mat'); %Load the Lie element generator

%Set the first template point:

%xhatp1 = [4.40984988833039e-4; 
		 %-4.39839482692374e-4; 
		  %0.728766350256459; 
		 %-0.130841217305973];
		 
[x, fval, info] = fsolve(@fxy, randn(4,1));

xhatp1 = x;

%Template tangent:
tp1 = T*xhatp1;

%Chart Border function
xhatstar = @(xhatstar1, xhatstar2, xhatp) [xhatstar1; 
										   xhatstar2; 
										   -(xhatstar2*(xhatp(2)*xhatp(3) - 2*xhatp(1)*xhatp(4)))/(4*(xhatp(3)^2 + xhatp(4)^2)) - (xhatstar1*(xhatp(1)*xhatp(3) + 2*xhatp(2)*xhatp(4)))/(4*(xhatp(3)^2 + xhatp(4)^2)); 
										   -(xhatstar1*(-2*xhatp(2)*xhatp(3) + xhatp(1)*xhatp(4)))/(4*(xhatp(3)^2 + xhatp(4)^2)) - (xhatstar2*(2*xhatp(1)*xhatp(3) + xhatp(2)*xhatp(4)))/(4*(xhatp(3)^2 + xhatp(4)^2))];

initial=-2;
final=2;
delta=0.1;
j = 1;
xhs=[1;1;1;1];
	
for xhatstar1 = initial:delta:final
	for xhatstar2 = initial:delta:final
			
				
				%Compute points on the chart borders:
				
				xhs(:,j) = xhatstar(xhatstar1, xhatstar2, xhatp1);
								
				j = j+1;
				%pause;
				
	end
end


plot3(xhs(1,:), xhs(2,:), xhs(3,:), '.r')
