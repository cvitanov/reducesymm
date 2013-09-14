clear all
hold off
clc

findzero = 1;

load xi.mat;
x(:,1) = xi;
i=1;
zero = 0;

xzero = [0.00154058117207164;
		-0.43983700573338302;
		 0.72941144405670832;
		-0.12719576522549556];

deltaphi = 0.01;
tol=4e-4;
phi=0;

if findzero
	
	while abs(xzero(1)) > tol
		
		phi = phi + deltaphi;
		xzero = LieElement(phi,xzero)
		
		if xzero(1) < 0
			
			xzero = LieElement(-phi,xzero)
			phi = phi - deltaphi;
			deltaphi = deltaphi / 2
	
			pause;
			
		end
		
	end

xzero

else

	for phi = 0.01:0.01:2*pi
		
		i = i+1;
		x(:,i) = LieElement(phi,xi);
		
		if x(1,i)*x(1,i-1)<0
		
			zero(length(zero)+1)=i;
		
		end
	
	end
	plot3(x(1,:), x(2,:), x(4,:));

end
