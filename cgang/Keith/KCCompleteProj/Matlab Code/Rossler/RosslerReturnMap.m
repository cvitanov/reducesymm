%Written by K. Carroll 5/1/2012
a = 0.2;
b = 0.2;
c = 5.7;
top = (c-sqrt(c^2-4*a*b))/2;
x0n = [top, -(top/a),(top/a)];
n = [1,0,0];
th = 2*pi - atan2(n(1:1,1), n(1:1,2))
numofstpos = 10;
numofcrossover = 100;
%number of starting positions on the unstable equilibrium
numofit = 20000;
% define two vectors one describing current pos, the other describing the
% previous
x1 = zeros(numofstpos,3);
x2 = zeros(numofstpos,3);
xtrack = zeros(numofstpos*numofcrossover, 5); 
sv = startingdirection(x0n);
for k=1:numofit
    x1(k:k,1:3) = x0n(1:1,1:3)+k*10^-6*(sv(1:3,1)');
end

for j=1:numofcrossover
    for i=1:numofstpos
        x2(i:i,1:3) = RKutta(@Rossler, x1(i:i,1:3),.01);
        while ~(crossover(x1(i:i,1:3), x2(i:i,1:3), n(1:1,1:3), x0n(1:1,1:3))) 
            x1(i:i,1:3) = x2(i:i,1:3);
            x2(i:i,1:3) = RKutta(@Rossler, x1(i:i,1:3),.01);
        end
        vadd = crossoverpos(x1(i:i,1:3), x2(i:i,1:3), n(1:1,1:3), x0n(1:1,1:3));
        xtrack(numofstpos*(j-1)+i:numofstpos*(j-1)+i,1:5) = [vadd, i, j];
        x1(i:i,1:3) = x2(i:i,1:3);
    end
    j
end

numofit = 5000;
x3 = zeros(numofit,3);
x3(1:1,1:3) = x0n(1:1,1:3)+(sv(1:3,1)')*10^-6;

for i=2:numofit
    x3(i:i,1:3) = RKutta(@Rossler, x3(i-1:i-1,1:3),.01);
end

dlmwrite('UnstableManifold.txt', xtrack, 'precision', 15);

hold all
plot3(x3(:,1), x3(:,2), x3(:,3), 'b');
plot3(xtrack(:,1), xtrack(:,2), xtrack(:,3), 'r', 'linestyle', 'none', 'Marker', 'o');
plot3(x0n(:,1), x0n(:,2), x0n(:,3), 'k', 'linestyle', 'none', 'Marker', 'o');
