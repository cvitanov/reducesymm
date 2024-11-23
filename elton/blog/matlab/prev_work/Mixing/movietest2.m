% trajectory_matlabRK4.m
% John Elton, updated 4/27/2008

% Create and plot a trajectory for a specified velocity field.
% Load the spectral coefficients and the geometry settings, then
% pick an initial point and a time length. 
% Calls 'u_sum' to create the physical velocity field, then uses
% Matlab's ode45 RK4 integration.

clear;

 global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N uspec1 uspec2 uspec3

%%%%%%%%%%%%%%   Setup    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

open('UB.mat');
uu = ans.UB;
open('UB_geom.mat');
ug = ans.UB_geom;

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
 
 t0 = 0;               % starting time
 tf = 40;              % final time
 Py = 4;          % number of trajectories (+1) to compute along y
 Pz = 4;          % number of trajectories (+1) to compute along z 
 
 tic
 count = 1;
 for jy = 1:Py-1
     for jz = 1:Pz-1
         r0 = [0 -1+(2*jy/Py) jz*(Lz/Pz)];   % initial position of jth particle        
         [t,r] = ode45('u_sum',[t0 tf],r0);   
         R{1,count} = r;
         if count == 1
            sizemin = size(r,1);
         elseif size(r,1) < sizemin
            sizemin = size(r,1);
         end
         count = count+1;
     end
 end
 toc
 
 for kk = 1:sizemin
    v1 = []; v2 = []; v3 = [];
    for j = 1:count-1
    v1 = [v1; R{1,j}(kk,1)];
    v2 = [v2; R{1,j}(kk,2)];
    v3 = [v3; R{1,j}(kk,3)];
    end
    plot3(v1,v3,v2,'LineStyle','-','Color',[0 0 1],'Marker','.')
    % mesh(v1,v3,v2)
   % plot3(R{1,:}(kk,1),R{1,:}(kk,3),R{1,:}(kk,2),'LineStyle','-','Color',[0 0 1],'Marker','.')
    grid on
    axis([0 Lx 0 Lz a b])
    xlabel('x = [0,Lx]'); zlabel('y = [-1,1]'); ylabel('z = [0,Lz]') 
    
    
    F(kk) = getframe;
 end

 movie(F,1,12)
 
 
 