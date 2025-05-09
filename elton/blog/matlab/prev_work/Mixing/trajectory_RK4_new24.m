% trajectory_matlabRK4.m
% John Elton, updated 12/11/24

% Calls 'u_sum' to create the physical velocity field, then uses
% Matlab's ode45 RK4 integration.


clear
global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N uspec1 uspec2 uspec3 nx ny nz dx dy dz v_grid r0 jy jz xcord rr

open('prev_work/Mixing/UB.mat');     % Upper Branch equilibrium
uu = ans.UB;
open('UB_geom.mat');
ug = ans.UB_geom;
open('v_g.mat');
v_grid = ans.v_g;


% geometry settings, don't edit
Mx = ug(1); My = ug(2); Mz = ug(3); Nx = ug(4); Ny = ug(5); Nz = ug(6);
Nd = ug(7); Lx = ug(8); Lz = ug(9); a = ug(10); b = ug(11);
lx = ug(12); lz = ug(13); alpha = ug(14); gamma = ug(15);
N = size(uu,1);
% create a column vector containing the spectral coefficients, for each component
uspec1 = zeros(N/3,1);
uspec1 = uu(1:N/3,1) + i*uu(1:N/3,2);
uspec2 = zeros(N/3,1);
uspec2 = uu(N/3+1:2*N/3,1) + i*uu(N/3+1:2*N/3,2);
uspec3 = zeros(N/3,1);
uspec3 = uu(2*N/3+1:N,1) + i*uu(2*N/3+1:N,2);
% end geometry settings


% define stagnation points
sp1_x = Lx/2; sp1_y = 0.0; sp1_z = Lz/4;
sp2_x = Lx/2; sp2_y = 0.0; sp2_z = (3*Lz)/4;
sp3_x = 0.0; sp3_y = 0.0; sp3_z = Lz/4;
sp4_x = 0.0; sp4_y = 0.0; sp4_z = (3*Lz)/4;
sp_n1_x = 2.35105561774981; sp_n1_y = 0.42293662349708; sp_n1_z = 0.65200166068573;
sp_n2_x = 3.16051044117966; sp_n2_y = -0.42293662349708; sp_n2_z = 0.60463540075018;


% % interpolation grid size
% nx = 48;           
% ny = 35;
% nz = 48;
% % increments
% dx = Lx/nx;
% dy = (b-a)/ny;
% dz = Lz/nz;



%%%%% Input parameters for run

% point = [sp1_x,sp1_y,sp1_z];
% open('evecs_eq1.mat');
% evecs_eq1 = ans.evecs_eq1;
% method = 'eigen';   % can be 'eigen', 'sphere'
% dim = 1;
% evecs = evecs_eq1;
% delt = 1/10;
% forward_time = false;
% color = 'r';
% num_t_evec = 10;          % number of trajectories (+1) to compute along evec
% rmax = 0.005;
% t_0 = 0;
% t_f = 100;
% t_step = 0.1;

% point = [sp1_x,sp1_y,sp1_z];
% open('evecs_eq1.mat');
% evecs_eq1 = ans.evecs_eq1;
% method = 'sphere';   % can be 'eig', 'sphere'
% dim = 2;
% evecs = evecs_eq1;
% delt = 1/100;
% forward_time = true;
% color = 'g';
% num_t_evec = 10;          % number of trajectories (+1) to compute along evec
% Py = 10;
% Pz = 10;
% rmax = 0.005;
% t_0 = 0;
% t_f = 100;
% t_step = 0.1;

point = [sp2_x,sp2_y,sp2_z];
open('evecs_eq2.mat');
evecs_eq2 = ans.evecs_eq2;
method = 'sphere';   % can be 'eig', 'sphere'
dim = 2;
evecs = evecs_eq2;
delt = 1/100;
forward_time = true;
color = 'k';
num_t_evec = 10;          % number of trajectories (+1) to compute along evec
Py = 10;
Pz = 10;
rmax = 0.01;
t_0 = 0;
t_f = 200;
t_step = 0.1;



% coords of initial point
% x_0 = 0.8; y_0 = -0.2; z_0 = 1.8;
% x_0 = 1.1; y_0 = -0.3; z_0 = 1.9;
% x_0 = 0.0; y_0 = -1.0; z_0 = 0.0;
% x_0 = 2.76244969613140; y_0 = 0.0; z_0 = 0.62086497079296;
x_0 = point(1); y_0 = point(2); z_0 = point(3);

np = (t_f-t_0)/t_step;

if forward_time
    tspan = [t_0:t_step:t_f];         
else
    tspan = [t_f:-t_step:t_0];     
end

if isequal(method, 'eig')
    if dim == 1
        for sides = 1:3   % (1 and 2)
            if sides == 1
                dr = delt*evecs(:,1);
            else
                dr = -delt*evecs(:,1);
            end
            r0 = [x_0+dr(1) y_0+dr(2) z_0+dr(3)]; 
            [t,r] = ode45('u_sum',tspan,r0);
            plot3(r(:,1),r(:,3),r(:,2),'LineStyle','-','Marker','.','Color',color)
            grid on
            offset = 0.5;
            axis([0-offset Lx+offset 0-offset Lz+offset a-offset/2 b+offset/2])
            xlabel('x = [0,Lx]'); zlabel('y = [-1,1]'); ylabel('z = [0,Lz]')      % y,z switch to match channelflow convention
            title('Upper Branch Equilibrium')
            plot3(x_0,z_0,y_0,'LineStyle','-','Marker','.','Color','k')  % marker for SP1
            hold on
        end
    elseif dim == 2
        for vec_num = (2:3)
            for jy = 1:num_t_evec-1
                dr = jy*(delt/10)*real(evecs(:,vec_num));
                r0 = [x_0+dr(1) y_0+dr(2) z_0+dr(3)];
                [t,r] = ode45('u_sum',tspan,r0);
                plot3(r(:,1),r(:,3),r(:,2),'LineStyle','-','Marker','.','Color',color)
                grid on
                offset = 0.5;
                axis([0-offset Lx+offset 0-offset Lz+offset a-offset/2 b+offset/2])
                xlabel('x = [0,Lx]'); zlabel('y = [-1,1]'); ylabel('z = [0,Lz]')      % y,z switch to match channelflow convention
                title('Upper Branch Equilibrium')
                plot3(x_0,z_0,y_0,'LineStyle','-','Marker','.','Color','k')  % marker for SP1
                hold on
            end
        end
    end

elseif isequal(method, 'sphere')           
    for jy = 1:Py-1
        xcord = rmax - 2*rmax*(jy-1)/(Py-2);
        rr = real(sqrt(rmax^2 - xcord^2));
        for jz = 1:Pz-1        
            r0 = [x_0+xcord y_0+rr*sin(2*pi*(jz-1)/(Py-2)) z_0+rr*cos(2*pi*(jz-1)/(Py-2))];    %%% parameterize with spherical coordinates           
            [t,r] = ode45('u_sum',tspan,r0);                    
            plot3(r(:,1),r(:,3),r(:,2),'LineStyle','-','Marker','.','Color',color)         
            grid on           
            offset = 0.5;
            axis([0-offset Lx+offset 0-offset Lz+offset a-offset/2 b+offset/2])
            xlabel('x = [0,Lx]'); zlabel('y = [-1,1]'); ylabel('z = [0,Lz]')      % y,z switch to match channelflow convention
            title('Upper Branch Equilibrium')
            plot3(x_0,z_0,y_0,'LineStyle','-','Marker','.','Color','k')  % marker for SP1
            hold on            
        end
    end
        
end




forward_time = true;
 
if forward_time
    tspan = [t_0:t_step:t_f];         
else
    tspan = [t_f:-t_step:t_0];     
end

 
tic
for jy = 1:Py-1
    xcord = rmax - 2*rmax*(jy-1)/(Py-2);
    rr = real(sqrt(rmax^2 - xcord^2));
    for jz = 1:Pz-1        
        r0 = [x_0+xcord y_0+rr*sin(2*pi*(jz-1)/(Py-2)) z_0+rr*cos(2*pi*(jz-1)/(Py-2))];    %%% parameterize with spherical coordinates
        %r0 = [x_0 y_0+(2*jy/Py) z_0+jz*(Lz/Pz)];   % initial position of jth particle   
        
        %[t,r] = ode45('u_interp',tspan,r0);
        [t,r] = ode45('u_sum',tspan,r0);
        % plot3(r0(1),r0(3),r0(2),'.')
        if forward_time
            plot3(r(:,1),r(:,3),r(:,2),'LineStyle','-','Marker','.','Color','g')
        else
            plot3(r(:,1),r(:,3),r(:,2),'LineStyle','-','Marker','.','Color','r')
        end 
        grid on
        % axis([0 Lx 0 Lz a b])
        offset = 0.5;
        axis([0-offset Lx+offset 0-offset Lz+offset a-offset/2 b+offset/2])
        xlabel('x = [0,Lx]'); zlabel('y = [-1,1]'); ylabel('z = [0,Lz]')      % y,z switch to match channelflow convention
        title('Upper Branch Equilibrium')
        hold on
        
    end
end
toc
 

 
 
 
 
 
 
 