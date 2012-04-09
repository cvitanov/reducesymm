% Calculates the group orbit X of a point x in the 
% 5D complex Lorenz system
function X = GroupOrbit(x)

theta = linspace(0,2*pi,50);

for j = 1:length(theta)
    g = gCLE(theta(j));  
    gx = g*x(:);
    X(j,1) = gx(1);
    X(j,2) = gx(2);
    X(j,3) = gx(3);
    X(j,4) = gx(4);
    X(j,5) = gx(5);
end
    