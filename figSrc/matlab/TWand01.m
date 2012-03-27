% Fig01.m by Daniel Borrero 03/23/2012
% Draws the group orbit of the TW_1 relative equilibrium and the relative
% periodic orbit 01 in the full space (no symmetry reduction)

clear all

% Part I: Calculate and plot group orbit of TW_1
% ------------------------------------------------------------------

% Coordinates of relative equilibrium TW_1
TWEQ = [8.484917778684904, -0.077135616169863, 8.485619011559175, 0, 26.999917355371903]; 
% Calculate group orbit of TW_1 by just rotating it by a series of finite
% angles between 0 and 2*pi
theta = linspace(0,2*pi,1000);

for j = 1:length(theta)
    g = zeros(5,5);
    g(1,1) = cos(theta(j));
    g(1,2) = sin(theta(j));
    g(2,1) = -sin(theta(j));
    g(2,2) = cos(theta(j));   
    g(3,3) = cos(theta(j));
    g(3,4) = sin(theta(j));
    g(4,3) = -sin(theta(j));
    g(4,4) = cos(theta(j));
    g(5,5) = 1; 
    
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
Tp = 1.542038909219227;

%Define initial condition for 01 relative periodic orbit
Xi = [2.38121811826785,14.1647812932010,2.26908302258214,14.1630085502492,34.9162551935267]'; 

for j = 1:length(theta)
    g = zeros(5,5);
    g(1,1) = cos(theta(j));
    g(1,2) = sin(theta(j));
    g(2,1) = -sin(theta(j));
    g(2,2) = cos(theta(j));   
    g(3,3) = cos(theta(j));
    g(3,4) = sin(theta(j));
    g(4,3) = -sin(theta(j));
    g(4,4) = cos(theta(j));
    g(5,5) = 1; 
    
    gy = g*Xi(:);
    
    Y(j,1) = gy(1);
    Y(j,2) = gy(2);
    Y(j,3) = gy(3);
    Y(j,4) = gy(4);
    Y(j,5) = gy(5);
end

% Define 5D complex Lorenz system
deq = inline('[-10*p(1)+10*p(3);-10*p(2)+10*p(4);(-p(5)+28).*p(1)-p(3)-0.1*p(4);(-p(5)+28).*p(2)+0.1*p(3)-p(4);-8/3*p(5)+p(1).*p(3)+p(2).*p(4)]','t','p');

% Integrate equations of motion from t = 0 to t = Tp
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
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
hold on
plot3(Y(:,1),Y(:,2),Y(:,5),'Color',[0 1 0],'LineWidth',2)

grid on

% Label axes
fsize = 14;
xlabel('x_1','Fontsize',fsize)
ylabel('x_2','Fontsize',fsize)
zlabel('z','Fontsize',fsize)

% Label orbits
legend('TW_1','01','Group Orbit of 01(t=0)')

% % Calculate shortest distance between group orbit and 01
% N = 100000000000;
% kmin = 1e6;
% jmin = 1e6;
% for j = 1:length(theta)
%     for k = 1:length(t)
%         if norm(X(j,:)-pp(k,:)) < N
%             N = norm(X(j,:)-pp(k,:));
%             kmin = k;
%             jmin = j;
%         end
%     end
% end
%         
% plot3(pp(kmin,1),pp(kmin,2),pp(kmin,5),'bo')
% plot3(X(jmin,1),X(jmin,2),X(jmin,5),'ro')
% 
% orbmin = pp(kmin,:)
% tmin = X(jmin,:)
