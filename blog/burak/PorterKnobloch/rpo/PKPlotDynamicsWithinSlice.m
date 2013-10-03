clear
hold off
clc

%Tutorial for making high quality octave plots:
%http://cresspahl.blogspot.com/2012/05/high-quality-plots-in-octave.html
%.eps outputs are converted to .pdf by epspdf

load dynamicswithintheslice.mat;
xhatslicedyn = xhat;
clear xhat;
load movingframes.mat;
xhatsymmred = xhat;
clear xhat;
load timeev.mat

figure(1)
subplot(1,2,1)
plot3(xhatslicedyn(1,:),xhatslicedyn(3,:),xhatslicedyn(4,:))
view(120,15)

xlabel('x_1')
ylabel('y_1')
zlabel('y_2')
title('within the slice')
view(120,15)

subplot(1,2,2)
plot3(xhatsymmred(1,:),xhatsymmred(3,:),xhatsymmred(4,:))
view(120,15)

xlabel('x_1')
ylabel('y_1')
zlabel('y_2')
title('on the moving frame')
view(120,15)

tarray = 1:size(phi,2);
tarray = tarray*deltat;

print('reddyn.eps','-S900,450',
'-depsc2',
'-F:Helvetica:9',
'-tight'
)

figure(2)

subplot(3,1,1)
plot(tarray,phi)
ylabel('\phi (t)')

tauarray = 1:size(phitau,2);
tauarray = tauarray*deltat;

subplot(3,1,2)
plot(tauarray, phitau)
ylabel('\phi ( \tau )')

subplot(3,1,3)

xdiff = xreconstructed - x;
plot(tarray, xdiff(1,:).^2+xdiff(2,:).^2+xdiff(3,:).^2+xdiff(4,:).^2)
xlabel('Time')
ylabel('|x_r (t) - x(t)|^2')

print('difference.eps','-S900,900',
'-depsc2',
'-F:Helvetica:9',
'-tight'
)

figure(3)

subplot(1,2,1)
plot3(xreconstructed(1,:), xreconstructed(3,:), xreconstructed(4,:));
xlabel('x_1')
ylabel('y_1')
zlabel('y_2')
view(120,30)
title('Reconstructed from the slice')

subplot(1,2,2)
plot3(x(1,:), x(3,:), x(4,:));
xlabel('x_1')
ylabel('y_1')
zlabel('y_2')
view(120,30)
title('Integrated in the full state space')

print('fullssp.eps','-S900,450',
'-depsc2',
'-F:Helvetica:9',
'-tight'
)
