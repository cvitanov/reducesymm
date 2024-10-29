clear;
P1 = 100;
P2 = 100;
R = 1;

function r = sphere_param(P1,P2,R)

for j1 = 1:P1
    for j2 = 1:P2
        x = R*sin(2*pi*(j1-1)/(P1-1))*cos(2*pi*(j2-1)/(P2-1));
        y = R*sin(2*pi*(j1-1)/(P1-1))*sin(2*pi*(j2-1)/(P2-1));
        z = R*cos(2*pi*(j1-1)/(P1-1));
        plot3(x,y,z,'.')
        grid on
        hold on
    end
end

hold off