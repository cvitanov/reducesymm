% trajectory_matlabRK4.m
% John Elton, updated 4/27/2008

% Calls 'u_sum' to create the physical velocity field, then uses
% Matlab's ode45 RK4 integration.


clear
global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N uspec1 uspec2 uspec3 nx ny nz dx dy dz v_grid r0 jy jz xcord rr

open('prev_work/Mixing/UB.mat');     % Upper Branch equilibrium
uu = ans.UB;
open('UB_geom.mat');
ug = ans.UB_geom;
tic
open('v_g.mat');
toc
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
 
nx = 48;           % interpolation grid size
ny = 35;
nz = 48;

dx = Lx/nx;
dy = (b-a)/ny;
dz = Lz/nz;
 
 Py = 8;          % number of trajectories (+1) to compute along y
 Pz = 8;          % number of trajectories (+1) to compute along z
 rmax = 0.1;
 
 for iii = 2:2
     if iii == 1 
     t0 = 110;               % starting time
     tf = 0;              % final time
     tspan = [t0:-0.5:tf];     
     else
     t0 = 0;               % starting time
     tf = 500;              % final time
     tspan = [t0:0.5:tf];   
     end
 
 tic
 for jy = 1:Py-1
     xcord = rmax - 2*rmax*(jy-1)/(Py-2);
     rr = real(sqrt(rmax^2 - xcord^2));
     for jz = 1:Pz-1        
        r0 = [0.8+xcord -.2+rr*sin(2*pi*(jz-1)/(Py-2)) 1.8+rr*cos(2*pi*(jz-1)/(Py-2))];    %%% parameterize with spherical coordinates
        %r0 = [xcord+1.1 -0.3+rr*sin(2*pi*(jz-1)/(Py-2)) 1.9+rr*cos(2*pi*(jz-1)/(Py-2))];
        %r0 = [0 -1+(2*jy/Py) jz*(Lz/Pz)];   % initial position of jth particle   
        %r0 = [2.76244969613140   0.00000000000000   0.62086497079296];
        

        %[t,r] = ode45('u_interp',tspan,r0);
        [t,r] = ode45('u_sum',tspan,r0);
       % plot3(r0(1),r0(3),r0(2),'.')
       if iii == 1
             plot3(r(:,1),r(:,3),r(:,2),'LineStyle','-','Marker','.','Color','r')
        else
             plot3(r(:,1),r(:,3),r(:,2),'LineStyle','-','Marker','.','Color','g')
        end 
        grid on
        axis([0 Lx 0 Lz a b])
        xlabel('x = [0,Lx]'); zlabel('y = [-1,1]'); ylabel('z = [0,Lz]')      % y,z switch to match channelflow convention
        title('Upper Branch Equilibrium')
        hold on
        
       end
   end
 toc
 
 end
 
 
 
 
 
 
 