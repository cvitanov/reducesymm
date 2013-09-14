clear
hold off;
clc

load('movingframes.mat');

figure(1)

plot3(xhat(1,:),xhat(3,:),xhat(4,:))
view(120,15)

xlabel('\hat{x}_1')
ylabel('\hat{y}_1')
zlabel('\hat{y}_2')
title('symmetry reduced')

%subplot(2,2,1)

%plot3(xhat(1,:),xhat(2,:),xhat(3,:))
%view(120,15)

%xlabel('\hat{x}_1')
%ylabel('\hat{x}_2')
%zlabel('\hat{y}_1')

%subplot(2,2,2)

%plot3(xhat(1,:),xhat(3,:),xhat(4,:))
%view(120,15)

%xlabel('\hat{x}_1')
%ylabel('\hat{y}_1')
%zlabel('\hat{y}_2')

%subplot(2,2,3)

%plot3(xhat(3,:),xhat(4,:),xhat(2,:))
%view(120,15)

%xlabel('\hat{y}_1')
%ylabel('\hat{y}_2')
%zlabel('\hat{x}_2')

%subplot(2,2,4)

%plot3(xhat(4,:),xhat(2,:),xhat(1,:))
%view(120,15)

%xlabel('\hat{y}_2')
%ylabel('\hat{x}_2')
%zlabel('\hat{x}_1')
