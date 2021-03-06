clear all;
clc;

load('generator.mat'); %Load the Lie element generator

%Set the first template point:

%xhatp1 = [4.40984988833039e-4; 
		 %-4.39839482692374e-4; 
		  %0.728766350256459; 
		 %-0.130841217305973];
		 
%[x, fval, info] = fsolve(@fxy, randn(4,1));

%xhatp1 = x;
xhatp1 = randn(4,1);

%Template tangent:
tp1 = T*xhatp1;

%Chart Border function
xhatstar = @(xhatstar1, xhatstar2, xhatp) [xhatstar1; 
										   xhatstar2; 
										   -(xhatstar2*(xhatp(2)*xhatp(3) - 2*xhatp(1)*xhatp(4)))/(4*(xhatp(3)^2 + xhatp(4)^2)) - (xhatstar1*(xhatp(1)*xhatp(3) + 2*xhatp(2)*xhatp(4)))/(4*(xhatp(3)^2 + xhatp(4)^2)); 
										   -(xhatstar1*(-2*xhatp(2)*xhatp(3) + xhatp(1)*xhatp(4)))/(4*(xhatp(3)^2 + xhatp(4)^2)) - (xhatstar2*(2*xhatp(1)*xhatp(3) + xhatp(2)*xhatp(4)))/(4*(xhatp(3)^2 + xhatp(4)^2))];

%Good choice inequalities:
%If both rows of the following vector is lesser than zero for the range given for chart borders, then the choce of templates is a good one.
 
%ineq = @(xhatp1, xhatp2, tp1, tp2, xhatstar1, xhatstar2) [dot(xhatp1-xhatp2, tp2) * dot(xhatstar1 - xhatp2, tp2);
%														  dot(xhatp2-xhatp1, tp1) * dot(xhatstar2 - xhatp1, tp1)];



%Transform the first template, pick that as the second template and check if it is a good one:

deltaphi = 0.01;
phi = -pi;
initial=-2;
final=2;
delta=1;
found=0;

while not(found)
	
	phi = phi+deltaphi;
	
	%xhatp2 = LieElement(phi, xhatp1); %Pick a second template candidate
	xhatp2 = randn(4,1);
	tp2 = T*xhatp2; %Second template tangent
	
	found = 1;
	j = 0;
	
	for xhatstar11 = initial:delta:final
		for xhatstar12 = initial:delta:final
			for xhatstar21 = initial:delta:final
				for xhatstar22 = initial:delta:final
				
					j = j + 1;
					%Compute points on the chart borders:
					
					xhatstar1 = xhatstar(xhatstar11, xhatstar12, xhatp1);
					xhatstar2 = xhatstar(xhatstar21, xhatstar22, xhatp2);
					
					if j == 1
						dot1=dot(xhatstar1, tp2);
						dot2=dot(xhatstar2, tp1);
					else
						dot1new=dot(xhatstar1, tp2);
						dot2new=dot(xhatstar2, tp1);
						
						if dot1*dot1new<0 || dot2*dot2new<0
							found =0;
							break;
						else
							dot1 = dot1new;
							dot2 = dot2new;
						end
					end
					%%Check the inequalities
					
					%ineqs = ineq(xhatp1, xhatp2, tp1, tp2, xhatstar1, xhatstar2);
					
					%j = j+1;
					%%pause;
					
					%if ineqs(1) >= 0 || ineqs(2) >= 0
						%found = 0;
						%break;
					%end
					
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
	
	if phi > pi - deltaphi
		break;
	end

end

found = found
%phi = phi
%fprintf('phi = %0.3f pi \n', phi./pi)
xhatp1 = xhatp1
xhatp2 = xhatp2
save('xhat.mat', 'xhatp1', 'xhatp2')
