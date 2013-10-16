clear all
hold off

clc

%The Rotation matrix operation:
R = @(theta, x) [cos(theta) -sin(theta);
				 sin(theta) cos(theta)]*x;

%Template point:
xhatp=[12;8;0;0];

%Angle between the Poincare section and the x axis:
theta = atan(8/12);

%Normal:
nhat=[R(pi/2, xhatp(1:2));0;0];

ps = poincare(xhatp,nhat,'psect.mat',-1);

n = size(ps,1) - 1; %get the dimensions

V=[xhatp/norm(xhatp) [0;0;1;0] nhat [0;0;0;1]];

psectrelative = ps(1:n, :) - xhatp;
psectprojected = [V' * psectrelative(1:n,:);
				  ps(n+1, :)];

%psectprojected = psect;
%psectprojected(1:2,:) = R(-theta,psect(1:2,:));

%r-Return map:

ri = sqrt(ps(1,:).^2 + ps(2,:).^2); 
riplus1 = ri(:,2:size(ri,2));
ri = ri(:,1:size(ri,2)-1);

figure(2)
plot3(psectprojected(1,:),psectprojected(2,:),psectprojected(4,:), '.', 'markersize', 5)
title('Poincare section')

figure(3)

%returnmap('psect.mat')
plot(ri,riplus1,'.', 'markersize', 5)
hold on
plot(min(ri):0.001:max(ri), min(ri):0.001:max(ri), 'k')
title('r-Return map')
xlabel('r_i')
ylabel('r_{i+1}')
