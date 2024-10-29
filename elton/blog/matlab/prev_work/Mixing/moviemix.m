% moviemix.m
% John Elton, updated May 05, 2008


keep v_grid
global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N uspec1 uspec2 uspec3 nx ny nz dx dy dz 

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
 
 nx = 144;           % interpolation grid size
ny = 105;
nz = 144;

dx = Lx/nx;
dy = (b-a)/ny;
dz = Lz/nz;
 
 P = 41;               % number of trajectories to compute
 t0 = 0;               % starting time
 tf = 25;              % final time
 tspan = [t0:0.5:tf];
 rr = 0.1;             % radius of IC circle
 
 tic
 for j = 1:P-1
  %  r0 = [(Lx/2)-1.1+j/10 0.2 Lz/2];
   r0 = [(Lx/2)+rr*cos(2*pi*(j-1)/(P-2)) rr*sin(2*pi*(j-1)/(P-2)) Lz/4];
  % r0 = [(Lx/2) rr*sin(2*pi*(j-1)/(P-2)) rr*cos(2*pi*(j-1)/(P-2))+3*Lz/4];   % initial position of jth particle  
    [t,r] = ode45('u_sum',tspan,r0);            
    R{1,j} = r;
 end
 toc
 
 sizer = size(r,1);
 for kk = 1:sizer
    v1 = []; v2 = []; v3 = [];
    for j = 1:P-1
        v1 = [v1; R{1,j}(kk,1)];
        v2 = [v2; R{1,j}(kk,2)];
        v3 = [v3; R{1,j}(kk,3)];
    end
  %  plot3(v1,v3,v2,'LineStyle','-','Color',[0 0 1])
    plot3(v1,v3,v2,'.')
    grid on
    axis([0 Lx 0 Lz a b])
    xlabel('x = [0,Lx]'); zlabel('y = [-1,1]'); ylabel('z = [0,Lz]')    
    F(kk) = getframe;
 end

 % save('movie5','F')
  movie(F,1,6)
 
 
 