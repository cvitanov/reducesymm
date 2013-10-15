clear
clc

%This m-file projects the flow on the Moving - Frames onto the 3D Gram-
%Schmidt basis

%Read the moving frames data:

load movingframes.mat

%Cartesian basis for full state space:

ex1 = [1;0;0;0];
ex2 = [0;1;0;0];
ey1 = [0;0;1;0];
ey2 = [0;0;0;1];

ex1GS = LieElement(pi/4,ex1);
ex2GS = LieElement(pi/4,ex2);
ey1GS = LieElement(pi/4,ey1);
ey2GS = LieElement(pi/4,ey2);

V = [ex1GS ex2GS ey1GS ey2GS];

for i = 1:size(xhat,2)

	for j = 1:size(V,2)
	
	xhatGS(j,i) = (V(:,j)'*xhat(:,i));
	
	end
	
end

xhatGS = xhatGS(2:4, :);

plot3(xhatGS(1,:), xhatGS(2,:), xhatGS(3,:))

view(285,25)

save('gramschmidt.mat','xhatGS')
