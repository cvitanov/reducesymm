% moviemix.m
% John Elton, updated May 05, 2008


clear;
global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N uspec1 uspec2 uspec3 nx ny nz dx dy dz v_grid

open('UB.mat');
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
 
 P = 20;               % number of trajectories to compute
 t0 = 0;               % starting time
 tf = 50;              % final time
 tspan = [t0:1:tf];
 h = 1;                % step size
 np = (tf-t0)/h;
 rmax = 0.1;
 
 tic
 for jy = 1:P-1
     xcord = rmax - 2*rmax*(jy-1)/(P-2);
     rr = sqrt(rmax^2 - xcord^2);
     for jz = 1:P-1        
        r = [1.1+xcord .85+rr*sin(2*pi*(jz-1)/(P-2)) 1.8+rr*cos(2*pi*(jz-1)/(P-2))];
       
        for i = 1:np
   
        k1 = u_interp(0,r(i,:));
        k2 = u_interp(0,r(i,:)+0.5*h*k1');
        k3 = u_interp(0,r(i,:)+0.5*h*k2');
        k4 = u_interp(0,r(i,:)+h*k3');
        r(i+1,:) = r(i,:) + (h/6)*(k1' + 2*k2' + 2*k3' + k4');


        %[t,r] = ode45('u_sum',tspan,r0); 
        %[t,r] = ode45('u_interp',tspan,r0);
        end
        R{1,jy,jz} = r;
       
     end
 end
 toc
 
 sizer = size(r,1);
 for kk = 1:sizer
    v1 = []; v2 = []; v3 = [];
    for jy = 1:P-1
        for jz = 1:P-1
        v1 = [v1; R{1,jy,jz}(kk,1)];
        v2 = [v2; R{1,jy,jz}(kk,2)];
        v3 = [v3; R{1,jy,jz}(kk,3)];
        end
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
 
 
 