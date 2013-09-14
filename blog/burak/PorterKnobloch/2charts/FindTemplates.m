function [results] = FindTemplates(xhatp1)

load('generator.mat'); %Load the Lie element generator

%Set the first template point:

%xhatp1 = [4.40984988833039e-4; 
		 %-4.39839482692374e-4; 
		  %0.728766350256459; 
		 %-0.130841217305973];

%Template tangent:
tp1 = T*xhatp1;

%Chart Border function
xhatstar = @(xhatstar1, xhatstar2, xhatp) [xhatstar1; 
										   xhatstar2; 
										   -(xhatstar2*(xhatp(2)*xhatp(3) - 2*xhatp(1)*xhatp(4)))/(4*(xhatp(3)^2 + xhatp(4)^2)) - (xhatstar1*(xhatp(1)*xhatp(3) + 2*xhatp(2)*xhatp(4)))/(4*(xhatp(3)^2 + xhatp(4)^2)); 
										   -(xhatstar1*(-2*xhatp(2)*xhatp(3) + xhatp(1)*xhatp(4)))/(4*(xhatp(3)^2 + xhatp(4)^2)) - (xhatstar2*(2*xhatp(1)*xhatp(3) + xhatp(2)*xhatp(4)))/(4*(xhatp(3)^2 + xhatp(4)^2))];

%Good choice inequalities:
%If both rows of the following vector is lesser than zero for the range given for chart borders, then the choce of templates is a good one.
 
ineq = @(xhatp1, xhatp2, tp1, tp2, xhatstar1, xhatstar2) [dot(xhatp1-xhatp2, tp2) * dot(xhatstar1 - xhatp2, tp2);
														  dot(xhatp2-xhatp1, tp1) * dot(xhatstar2 - xhatp1, tp1)];

%Transform the first template, pick that as the second template and check if it is a good one:

deltaphi = 0.01;
phi = 0;
initial=-2;
final=2;
delta=1;
found=0;

while not(found)
	
	phi = phi+deltaphi;
	
	xhatp2 = LieElement(phi, xhatp1); %Pick a second template candidate
	tp2 = T*xhatp2; %Second template tangent
	
	found = 1;
	j = 0;
	
	for xhatstar11 = initial:delta:final
		for xhatstar12 = initial:delta:final
			for xhatstar21 = initial:delta:final
				for xhatstar22 = initial:delta:final
				
					
					%Compute points on the chart borders:
					
					xhatstar1 = xhatstar(xhatstar11, xhatstar12, xhatp1);
					xhatstar2 = xhatstar(xhatstar21, xhatstar22, xhatp2);
	
					%Check the inequalities
					
					ineqs = ineq(xhatp1, xhatp2, tp1, tp2, xhatstar1, xhatstar2);
					
					j = j+1;
					%pause;
					
					if ineqs(1) >= 0 || ineqs(2) >= 0
						found = 0;
						break;
					end
					
				end
				if not(found)
					break;
				end
			end
			if not(found)
				break;
			end
		end
		if not(found)
			break;
		end		
	end
	
	if phi > 2*pi - deltaphi
		break;
	end

end

phi = phi
fprintf('phi = %0.3f pi \n', phi./pi)
xhatp2 = xhatp2
