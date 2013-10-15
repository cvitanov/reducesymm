clear all
hold off
figure(2)
clc

%The Rotation matrix operation:
R = @(theta, x) [cos(theta) -sin(theta);
				 sin(theta) cos(theta)]*x;

%Template point:
xhatp=[3.5;14;0;0];

%Angle between the Poincare section and the x axis:
theta = atan(14/3.5);

%Normal:
nhat=[R(pi/2, xhatp(1:2));0;0];

psect = poincare(xhatp,nhat);

psectprojected = psect;
psectprojected(1:2,:) = R(-theta,psect(1:2,:));

%r-Return map:

ri = sqrt(psect(1,:).^2 + psect(2,:).^2); 
riplus1 = ri(:,2:size(ri,2));
ri = ri(:,1:size(ri,2)-1);

subplot(1,2,1)
plot(psectprojected(1,:),psectprojected(3,:), '.', 'markersize', 5)
title('Poincare section')

subplot(1,2,2)
plot(ri,riplus1,'.', 'markersize', 5)
title('r-Return map')
xlabel('r_i')
ylabel('r_{i+1}')
