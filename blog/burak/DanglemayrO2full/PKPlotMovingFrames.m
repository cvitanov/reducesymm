clear
hold off;
clc

load('movingframes.mat');
load('xinvsubs.mat');
figure(1)

plot3(xhat(1,:),xhat(2,:),xhat(4,:))
view(120,15)

xlabel('x_1')
ylabel('x_2')
zlabel('y_2')
title('Symmetry Reduced')

hold on;

plot3(xinvsubs(1,:), xinvsubs(2,:), xinvsubs(4,:), 'r.')

print -dpdf symmetryreduced.pdf

xhatGS = LieElement(pi/4,xhat);

figure(2)

plot3(xhatGS(1,:), xhatGS(3,:), xhatGS(4,:))
xlabel('x_{GS,1}')
ylabel('y_{GS,1}')
zlabel('y_{GS,2}')
title('On Gram-Schmidt Coordinates')

print -dpdf gramschmidt.pdf
