%Written by K. Carroll 5/1/2012
b = 8/3;
p1 = 28;
e = 1/10;
sig = 10;
d = 1+e^2/((sig+1)^2);
r = sqrt(b*(p1-d));
r2 = sqrt(d)*r;
z = p1-d;
theta = -acos(1/sqrt(d));
x0n = [r*cos(theta), r*sin(theta), r2, 0, z];
x0s = [1,0,0,0,0];
rdis = dot(x0n(1:1,1:4),x0n(1:1,1:4));
T =zeros(5,5);
T(1:1,2:2) = -1;
T(2:2,1:1) = 1;
T(3:3,4:4) = -1;
T(4:4,3:3) = 1;
tx1 = T*(x0s');

opp = x0n(1:1,1:5)*tx1;
adj = x0n(1:1,1:4)*x0s(1:1,1:4)';
phi = atan(opp/adj);
d = sqrt(opp^2+adj^2);
c = adj/d;
s = opp/d;
x0nr(1:1,1:5) = [c*x0n(1:1,1:1)+s*x0n(1:1,2:2),-s*x0n(1:1,1:1)+c*x0n(1:1,2:2), c*x0n(1:1,3:3)+s*x0n(1:1,4:4),  -s*x0n(1:1,3:3)+c*x0n(1:1,4:4), x0n(1:1,5:5)];

n = [1, 0, 0, 0, 0];
numofstpos = 50;
% %number of starting positions on the unstable equilibrium
numofit = 15000;
% % define two vectors one describing current pos, the other describing the
% % previous
x1 = zeros(numofstpos, 5); 
x2 = zeros(numofstpos, 5);
x1r = zeros(numofstpos, 5); 
x2r = zeros(numofstpos, 5);
numofcross = 300;
% Define a vector which has the pos1, pos2,...posn, trajectory, crossing
xcross = zeros(numofcross*numofstpos, 7);
UnstableMan = zeros(numofstpos,606);
col = zeros(numofstpos,1);
sv = startingdirectionCLE(x0n);

for k=1:numofstpos
    x1(k:k,1:5) = x0n(1:1,1:5)+k*10^-6*(sv(1:5,1)');
    opp = x1(1:1,1:5)*tx1;
    adj = x1(1:1,1:4)*x0s(1:1,1:4)';
    phi = atan(opp/adj);
    d = sqrt(opp^2+adj^2);
    c = adj/d;
    s = opp/d;
    x1r(k:k,1:5) = [c*x1(k:k,1:1)+s*x1(k:k,2:2),-s*x1(k:k,1:1)+c*x1(k:k,2:2), c*x1(k:k,3:3)+s*x1(k:k,4:4),  -s*x1(k:k,3:3)+c*x1(k:k,4:4), x1(k:k,5:5)];
end

for i=1:numofcross
    for k=1:numofstpos
        x2(k:k,1:5) = RKutta(@ComplexLorenz, x1(k:k,1:5),.01);
        opp = x2(k:k,1:5)*tx1;
        adj = x2(k:k,1:4)*x0s(1:1,1:4)';
        phi = atan(opp/adj);
        d = sqrt(opp^2+adj^2);
        c = adj/d;
        s = opp/d;
        x2r(k:k,1:5) = [c*x2(k:k,1:1)+s*x2(k:k,2:2),-s*x2(k:k,1:1)+c*x2(k:k,2:2), c*x2(k:k,3:3)+s*x2(k:k,4:4),  -s*x2(k:k,3:3)+c*x2(k:k,4:4), x2(k:k,5:5)];
        while ~(crossover(x1r(k:k,1:5),x2r(k:k,1:5), n(1:1,1:5), x0nr(1:1,1:5)))
            x1(k:k,1:5) = x2(k:k,1:5);
            x1r(k:k,1:5) = x2r(k:k,1:5);
            x2(k:k,1:5) = RKutta(@ComplexLorenz, x1(k:k,1:5),.01);
            opp = x2(k:k,1:5)*tx1;
            adj = x2(k:k,1:4)*x0s(1:1,1:4)';
            phi = atan(opp/adj);
            d = sqrt(opp^2+adj^2);
            c = adj/d;
            s = opp/d;
            x2r(k:k,1:5) = [c*x2(k:k,1:1)+s*x2(k:k,2:2),-s*x2(k:k,1:1)+c*x2(k:k,2:2), c*x2(k:k,3:3)+s*x2(k:k,4:4),  -s*x2(k:k,3:3)+c*x2(k:k,4:4), x2(k:k,5:5)];        
        end
        vadd = crossoverpos(x1r(k:k,1:5),x2r(k:k,1:5), n(1:1,1:5), x0nr(1:1,1:5));
        xcross((i-1)*numofstpos+k:(i-1)*numofstpos+k, 1:7) = [vadd, k, i];
        x1(k:k,1:5) = x2(k:k,1:5);
        x1r(k:k,1:5) = x2r(k:k,1:5);
    end
    i
end

x3 = zeros(numofit,5);
x3r = zeros(numofit,5);
x3(1:1,1:5) = x0n(1:1,1:5)+(sv(1:5,1)')*10^-6;
opp = -x3(1:1,1:5)*tx1;
adj = x3(1:1,1:4)*x0s(1:1,1:4)';
phi = atan(opp/adj);
d = sqrt(opp^2+adj^2);
c = adj/d;
s = opp/d;
x3r(1:1,1:5) = [c*x3(1:1,1:1)+s*x3(1:1,2:2),-s*x3(1:1,1:1)+c*x3(1:1,2:2), c*x3(1:1,3:3)+s*x3(1:1,4:4),  -s*x3(1:1,3:3)+c*x3(1:1,4:4), x3(1:1,5:5)];

for i=2:numofit
    x3(i:i,1:5) = RKutta(@ComplexLorenz, x3(i-1:i-1,1:5),.01);
    opp = x3(i:i,1:5)*tx1;
    adj = x3(i:i,1:4)*x0s(1:1,1:4)';
    phi = atan(opp/adj);
    d = sqrt(opp^2+adj^2);
    c = adj/d;
    s = opp/d;
    x3r(i:i,1:5) = [c*x3(i:i,1:1)+s*x3(i:i,2:2),-s*x3(i:i,1:1)+c*x3(i:i,2:2), c*x3(i:i,3:3)+s*x3(i:i,4:4),  -s*x3(i:i,3:3)+c*x3(i:i,4:4), x3(i:i,5:5)];
end
dlmwrite('UnstableManifold.txt', xcross, 'precision', 15);

hold all
plot3(x3r(:,1), x3r(:,2), x3r(:,3), 'b');
plot3(xcross(:,1), xcross(:,2), xcross(:,3), 'r', 'linestyle', 'none', 'Marker', 'o');
plot3(x0nr(:,1), x0nr(:,2), x0nr(:,3), 'k', 'linestyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
