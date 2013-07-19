% Written by K. Carroll 5/1/2012
load startingpoints.txt
load ReducedPO.txt

s1 = size(startingpoints);
s2 = size(ReducedPO())
RPO = ReducedPO()
xi = startingpoints();

a = 0.2;
b = 0.2;
c = 5.7;
top = (c-sqrt(c^2-4*a*b))/2;
top2 = (c+sqrt(c^2-4*a*b))/2;

x0n = [top, -(top/a),(top/a)];
x0f = [top2, -(top2/a),(top2/a)];
n = s2(1:1,2);
n2 = s2(1:1,1);
numofit = n*1000;
M = zeros(n2, 4+s2(1:1,2));
x = zeros(numofit, 3);
xpolar = zeros(numofit, 4);
for k=1:n2
    xpolar(1:1,1:4) = [sqrt(xi(k:k,1)^2+xi(k:k,2)^2), atan2(xi(k:k,2), xi(k:k,1)), xi(k:k,3), 0];

    x(1:1,1:3) = xi(k:k,1:3);
    Jac = eye(3);
    dtheta = 2*pi/1000;
    for i=2:numofit
        xpolar(i:i,1:4) = RKutta(@rthetaz, xpolar(i-1:i-1,1:4), dtheta);
        x(i:i,1:3) = [xpolar(i:i,1)*cos(xpolar(i:i,2)),xpolar(i:i,1)*sin(xpolar(i:i,2)), xpolar(i:i,3)];
        dt = xpolar(i:i,4)-xpolar(i-1:i-1,4);
        Jac = Jacob(Jac, x(i:i,1:3), dt);
    end
    L = eig(Jac);
    T = xpolar(numofit:numofit,4);
    M(k:k,:) = [L', T, RPO(k:k,:)];
end
hold all
plot3(x(:,1), x(:,2), x(:,3), 'b');
plot3(x0n(:,1), x0n(:,2), x0n(:,3), 'k', 'linestyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
% zlim([0, 25]);
% axis off
xlabel('x')
ylabel('y')
zlabel('z')
view(-70, 14);
dlmwrite('values7it.txt', M, 'precision', 6)
