% Fig01.m by Daniel Borrero 03/23/2012
% Draws the group orbit of the TW_0 relative equilibrium and the relative
% periodic orbit 01 in the full space (no symmetry reduction)

clear all

% Part I: Calculate and plot group orbit of TW_0
% ------------------------------------------------------------------

% Coordinates of relative equilibrium TW_0
TWEQ = [8.4849, 0.077135, 8.4856, 0, 26.9989]; 

% Calculate group orbit of TW_0 by just rotating it by a series of finite
% angles between 0 and 2*pi
theta = linspace(0,2*pi,1000);

for j = 1:length(theta)
    g = gCLE(theta(j));  
    gx = g*TWEQ(:);
    X(j,1) = gx(1);
    X(j,2) = gx(2);
    X(j,3) = gx(3);
    X(j,4) = gx(4);
    X(j,5) = gx(5);
end


% Part II: Calculate and plot RPO 01 in the full space
% -----------------------------------------------------------------

% Define period for 01 relative periodic orbit
Tp = 1.54203890921923; 

%Define initial condition for 01 relative periodic orbit
Xi = [2.38121811826785,14.1647812932010,2.26908302258214,14.1630085502492,34.9162551935267]'; 

% Define 5D complex Lorenz system
deq = inline('[-10*p(1)+10*p(3);-10*p(2)+10*p(4);(-p(5)+28).*p(1)-p(3)-0.1*p(4);(-p(5)+28).*p(2)+0.1*p(3)-p(4);-2.7*p(5)+p(1).*p(3)+p(2).*p(4)]','t','p');

% Integrate equations of motion from t = 0 to t = Tp
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,pp] = ode45(deq,[0 Tp],Xi,options);
 
% Part 3: Plot results
% ---------------------------------------------------------------
figure(1)

% Plots the first(x1), second(x2) and fifth(z) coordinates of each trajectory
% Modify pp(:,n) to get different coordinates. Make sure to adjust axis
% labels
plot3(X(:,1),X(:,2),X(:,5),'Color',[1 0 0],'LineWidth',2)
hold on
plot3(pp(:,1),pp(:,2),pp(:,5),'Color',[0 0 1],'LineWidth',2)

grid on

% Label axes
fsize = 14;
xlabel('x_1','Fontsize',fsize)
ylabel('x_2','Fontsize',fsize)
zlabel('z','Fontsize',fsize)

% Label orbits
legend('TW_0','01')


