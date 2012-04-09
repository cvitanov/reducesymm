% Wurst01.m by Daniel Borrero 4/9/2012
% Draws full time/group orbit for 01 relative periodic orbit of 
% the Complex Lorenz flow. See comments below for details.

clear all

figure()

% Set Complex Lorenz flow parameters
rho1 = 28;
rho2 = 0;
sigma =10;
b = 8/3;
e = 1/10;

params = {rho1,rho2,sigma,b,e};


TWEQ = [ 8.48490000000000;... % Relative equilibrium as per blog discussion 
        -0.0771349999999995;...  % on 2012-03-25/26
         8.48560000000000;...
         0;...
         26.9990000000000;]

% Move equilibrium around group orbit by -pi/2 (looks better on final plot)
g = gCLE(-pi/2);
TWEQ = g*TWEQ;

Tp = 1.54203890921923; % period of 01 orbit from Siminos
THETA = [2.95261721254045;]; % phase shift of 01 orbit from Siminos
Xi = [2.38121811826785,...  % initial condition for 01 orbit from Siminos
      14.1647812932010,...
      2.26908302258214,...
      14.1630085502492,...
      34.9162551935267]';

% Generate and plot group orbits 
GX = GroupOrbit(Xi); % group orbit of initial condition for 01
GTW = GroupOrbit(TWEQ); % group orbit of relative equilibrium TW_1

% Plot initial conditions and it's group orbit for 01 orbit
% (set colors with [r g b] values where r, g, b are between 0 and 1
% in MarkerFaceColor and Color inputs)
plot3(GX(:,1),GX(:,2),GX(:,5),'--b','LineWidth',2),hold on
plot3(Xi(1),Xi(2),Xi(5),'Marker','o','MarkerFaceColor',[0 0 1], 'Color', [0 0 1],'MarkerSize',8)

% Plot TW_1 and its group orbit
plot3(GTW(:,1),GTW(:,2),GTW(:,5),'--r','LineWidth',2)
plot3(TWEQ(1),TWEQ(2),TWEQ(5),'Marker','o','MarkerFaceColor',[1 0 0], 'Color', [1 0 0],'MarkerSize',8)

% Calculate N = n+1 (16) periods of the 01 relative periodic orbit,
% making sure to avoid exponential blow up advancing the phase by
% THETA once per period.
for n = 0:15
    % Setup initial condition with phase advance
    g = gCLE(n*THETA);
    xi = g*Xi;
    % Set integrator options
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    % Integrate equations of motion
    [T,Y] = ode45(@ComplexLorenzEOM,[0 Tp],xi,options,params);
    
    % For first period, plot orbit with solid line
    if n == 0
        plot3(Y(:,1),Y(:,2),Y(:,5),'Color',[0 0 1],'Linewidth',2)
        plot3(Y(end,1),Y(end,2),Y(end,5),'Marker','o','MarkerFaceColor',[0 0 1], 'Color', [0 0 1],'MarkerSize',8)
    % For other periods, plot faint dotted lines
    else
        plot3(Y(:,1),Y(:,2),Y(:,5),'Color',[0 0.5 0.5],'Linewidth',1,'LineStyle',':')
    end
    drawnow
end


% Calculate and plot group orbit for a few select points on 01
% (may be confusing, but helps enhance 3D effect of surface drawn later)
% Comment this section to skip this step.
for j = 1:250:length(Y(:,1))
   GX = GroupOrbit(Y(j,1:5));
   plot3(GX(:,1),GX(:,2),GX(:,5),'Color',[0 .5 .5],'Linewidth',1,'LineStyle',':')
end


% Calculate distances from z-axis for select point on 01 orbit
% and extract their z-positions as a vector
R = sqrt(Y(1:150:end,1).^2 + Y(1:150:end,2).^2);
Z = Y(1:150:end,5);

% Caculate the surface coordinates for cylinders having the radii
% saved in R
[a,b,c] = cylinder(R);

% Generate matrix of z positions for points calculated on the cylinders
% calculated above
C = zeros(size(a));
for i = 1:length(a(:,1))
   C(i,:) = Z(i);
end

% Plot group orbit for entire group orbit of 01 relative periodic orbit.
surf(a,b,C,'FaceColor',[0 .9 .9],'EdgeColor','none','Linewidth',0.5,'LineStyle',':')
alpha(0.12) % Set transparency value for surface (values between 0.1 and 0.2 seem to work well)

% Set camera view. This can be adjusted by hand once the figure is rendered 
% but lighting may have to be adjusted for it to look good.
rotate3d
view([115 12]) % [Azymuthal rotation, elevation angle]
drawnow

% Set axes limits and hide tick marks
axis([-20 20 -20 20 0 60])
axis off

% Draw 3D arrows for axes. Nicer than regular quiver3 arrows
% Can be commented out if axes are going to be drawn with outside software
mArrow3([0 0 0], [0 0 55]) % Draw z axis
mArrow3([0 0 0], [20 0 0]) % Draw x1 axis
mArrow3([0 0 0], [0 20 0]) % Draw x2 axis

% Draw guidelines to draw axes with outside software 
% plot3([0 0],[0 0], [0 55],'k-','LineWidth',0.5) % Draw z axis guideline 
% plot3([0 20],[0 0], [0 0],'k-','LineWidth',0.5) % Draw x1 axis guideline 
% plot3([0 0],[0 20], [0 0],'k-','LineWidth',0.5) % Draw x2 axis guideline 
% drawnow

% Set lighting options (there may be more playing to do here
camlight headlight % Haven't played with this much. For more complicated lighting see next line
% camlight(115,12,'infinite'); % Set location of lighting source (see help camlight for more options) 
lighting phong % Set rendering style (see help lighting)
drawnow 
camlight headlight % Not sure why but makes lighting a little brighter
lighting phong