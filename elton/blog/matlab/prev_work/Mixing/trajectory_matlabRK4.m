% trajectory_matlabRK4.m
% John Elton, updated 5/23/2008

% Calls 'u_sum' to create the physical velocity field, then uses
% Matlab's ode45 RK4 integration.

keep v_grid evecs_eq2 sp5 sp6

global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N uspec1 uspec2 uspec3 v_grid dx dy dz  

%global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N uspec1 uspec2 uspec3 nx ny nz dx dy dz v_grid r0 jy jz xcord rr


open('UB.mat');      % Upper Branch equilibrium
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
 
nx = 5;           % interpolation grid size, one less than dimension of array
ny = 4;
nz = 5;

dx = Lx/nx;
dy = (b-a)/ny;
dz = Lz/nz;

 for iii = 2:2
     if iii == 2 
     t0 = 190;               % starting time
     tf = 0;              % final time
     tspan = [t0:-0.5:tf];     
     else
     t0 = 0;               % starting time
     tf = 150;              % final time
     tspan = [t0:0.5:tf];   
     end
 
        P1 = 41;         
        P2 = 8;          
        R = 0.01;
      
%         for j1 = 1:P1
%             for j2 = 1:P2
                  % sphere of ICs  
%                 x = R*sin(2*pi*(j1-1)/(P1-1))*cos(2*pi*(j2-1)/(P2-1));
%                 y = R*sin(2*pi*(j1-1)/(P1-1))*sin(2*pi*(j2-1)/(P2-1));
%                 z = R*cos(2*pi*(j1-1)/(P1-1));
%                 r0 = [x+Lx/2 y z+3*Lz/4];

%                   %  complex pair   eq2 eq3 eq5 eq6 
%                   v1 = real(evecs_eq5(:,2));
%                   v1 = (1/norm(v1))*v1;
%                   v2 = imag(evecs_eq5(:,2));
%                   v2twidle = v2 - dot(v2,v1)*v1;
%                   v2 = (1/norm(v2twidle))*v2twidle;
                  
%                   % unstalbe, real  eq1 or eq4
%                   v1 = evecs_eq1(:,2);
%                   v1 = (1/norm(v1))*v1;
%                   v2 = evecs_eq1(:,3);
%                   v2twidle = v2 - dot(v2,v1)*v1;
%                   v2 = (1/norm(v2twidle))*v2twidle;
                  
                    eps = 10^2;
                    k1 = 1;

%                   for j1 = 1:P1
%                       
%                       % set initial conditions along evecs
%                       if j1 < P1/2
%                       r0 = [Lx/2+(j1/eps)*v1(1) (j1/eps)*v1(2) 3*Lz/4+(j1/eps)*v1(3)];
%                       else
%                        r0 = [Lx/2+(k1/eps)*v2(1) (k1/eps)*v2(2) 3*Lz/4+(k1/eps)*v2(3)];
%                        k1 = k1+1;
%                       end

%                   % set initial conditions in a circle
%                   r0 = cos(2*pi*(j1-1)/(P1-1))*v1 + sin(2*pi*(j1-1)/(P1-1))*v2;
%                   r0 = [(1/100)*r0(1)+sp5(1) (1/100)*r0(2)+sp5(2) (1/100)*r0(3)+sp5(3)];

                     % for 1D manifolds
                   r0 = -(1/50)*evecs_eq2(:,1);
                   r0 = [r0(1)+Lx/2 r0(2) r0(3)+3*Lz/4]; 
%                    r0 = [r0(1) r0(2) r0(3)+Lz/4];
                   
                  tic  
                  %[t,r] = ode45('u_interp',tspan,r0);
                  [t,r] = ode45('u_sum',tspan,r0);
                  toc
                  
                  %plot3(r0(1),r0(3),r0(2),'.')
                  
                  if iii == 1
                  %plot3(r(:,1),r(:,3),r(:,2),'LineStyle','-','Marker','.','Color','b')
                  plot3(r(:,1),r(:,3),r(:,2),'LineStyle','-','Color','b')
                  else
                  plot3(r(:,1),r(:,3),r(:,2),'LineStyle','-','Color','k')
                  end 
                  grid on
                  axis([0 3*Lx/2 0 Lz a b])
                  xlabel('x = [0,Lx]'); zlabel('y = [-1,1]'); ylabel('z = [0,Lz]')      % y,z switch to match channelflow convention
                  title('Upper Branch Equilibrium')
                  hold on 
                  
%                   end
 
% %             end
% %         end
%        
  
 
 end
 
 plot3(Lx/2,3*Lz/4,0,'LineStyle','-','Marker','.','Color','k')
 plot3(Lx/2,Lz/4,0,'LineStyle','-','Marker','.','Color','k')
%  plot3(Lx,3*Lz/4,0,'LineStyle','-','Marker','.','Color','k')
%  plot3(Lx,Lz/4,0,'LineStyle','-','Marker','.','Color','k')
 plot3(0,3*Lz/4,0,'LineStyle','-','Marker','.','Color','k')
 plot3(0,Lz/4,0,'LineStyle','-','Marker','.','Color','k')
 plot3(sp5(1),sp5(3),sp5(2),'LineStyle','-','Marker','.','Color','k')
 plot3(sp6(1),sp6(3),sp6(2),'LineStyle','-','Marker','.','Color','k')