function [xhat, phitau]=findxhat(x, xhatp)

%This function takes the state space point x and finds its corresponding
%point on the slice defined by xhatp

load('generator.mat'); %Load the Lie element generator

%Template tangent
tp = T * xhatp;

%Slice condition:
SliceCondition = @(phi) dot(tp, LieElement(-phi, x)); %eq 10.44

deltaphi = 0.1; %starting phi-step

for phi=0:deltaphi:2*pi
		
	if SliceCondition(phi) < 0 && SliceCondition(phi+deltaphi)>0
					
		[phii, fval, info] = fsolve(SliceCondition, phi);
		
		phitau = phii;
		xhat = [LieElement(-phii, x); phitau];
			
		break;
			
	end
	
end
