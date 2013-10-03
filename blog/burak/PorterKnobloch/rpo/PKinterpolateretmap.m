clear
clc

load('psect1.mat'); %Load the Poincare section and the return map data.

[snsorted,index]=sort(sn); %find the peak

snplus1sorted = snplus1(index);

imax = find(snplus1sorted == max(snplus1sorted));

%interpolate in two pieces:

step = 1e-6;

x1 = snsorted(1:imax);
y1 = snplus1sorted(1:imax);
xinterp1 = snsorted(1):step:snsorted(imax);
yinterp1 = interp1(x1, y1, xinterp1, 'spline');

x2 = snsorted(imax:length(snsorted));
y2 = snplus1sorted(imax:length(snsorted));
xinterp2 = snsorted(imax):step:snsorted(length(snsorted));
yinterp2 = interp1(x2, y2, xinterp2, 'spline');

%clf();
%hold on
%plot(snsorted,snplus1sorted, '.')
%plot(xinterp1, yinterp1, 'k')
%plot(xinterp2, yinterp2, 'k')
%plot(min(snplus1):step:max(snplus1),min(snplus1):step:max(snplus1),'r')

interp2 = [xinterp2;
		   yinterp2];

tol = 1e-6;

irpo=find(abs(interp2(1,:)-interp2(2,:))<tol);
srpo=yinterp2(irpo);


idummy = find(y2<srpo);
idummy1 = idummy(1);
idummy2 = idummy(1)-1;

i1snplus1 = find(snplus1==y2(idummy1));
i2snplus1 = find(snplus1==y2(idummy2));

i1ps = i1snplus1 + 1;
i2ps = i2snplus1 + 1;

%Two Poincare section points between which the relative equilibrium lies:
x1=ps(1:3,i1ps)
x2=ps(1:3,i2ps)

xrpoGS = x1 + ((srpo - snplus1(i1snplus1))/(snplus1(i2snplus1) - snplus1(i1snplus1)))*(x2-x1);
xrpo = GS2x(xrpoGS)
