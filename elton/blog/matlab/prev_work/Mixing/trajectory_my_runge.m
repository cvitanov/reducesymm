keep v_grid evecs_eq2
global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N uspec1 uspec2 uspec3 dx dy dz v_grid 

%%%%%%%%%%%%%%   Setup    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 open('UB.mat');
uu = ans.UB;
open('UB_geom.mat');
ug = ans.UB_geom;
 Py = 11;          % number of trajectories (+1) to compute along y
 Pz = 11;          % number of trajectories (+1) to compute along z
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
 
 ti = 100;               % starting time
 tf = 0;              % final time
 h = -.1;                % step size
 np = (tf-ti)/h;
 rmax = 0.001;
 
 nx = 143;
 ny = 104;
 nz = 143;
 
 dx = Lx/nx;
 dy = (b-a)/ny;
 dz = Lz/nz;
 
%  for jy = 1:Py-1
%      xcord = rmax - 2*rmax*(jy-1)/(Py-2);
%      rr = sqrt(rmax^2 - xcord^2);
%      for jz = 1:Pz-1        
%         r0 = [Lx/2+xcord rr*sin(2*pi*(jz-1)/(Py-2)) Lz/4+rr*cos(2*pi*(jz-1)/(Py-2))];
 
         r0 = -(1/10)*evecs_eq2(:,1);
                  r0 = [r0(1)+Lx/2 r0(2) r0(3)+3*Lz/4]; 
r = [r0];
tic
for i = 1:np
    
% k1 = u_interp_c(0,r(i,:),Lx,Lz,[nx ny nz],v_grid);
% k2 = u_interp_c(0,r(i,:)+0.5*h*k1',Lx,Lz,[nx ny nz],v_grid);
% k3 = u_interp_c(0,r(i,:)+0.5*h*k2',Lx,Lz,[nx ny nz],v_grid);
% k4 = u_interp_c(0,r(i,:)+h*k3',Lx,Lz,[nx ny nz],v_grid);
    
   
% k1 = u_interp(0,r(i,:));
% k2 = u_interp(0,r(i,:)+0.5*h*k1');
% k3 = u_interp(0,r(i,:)+0.5*h*k2');
% k4 = u_interp(0,r(i,:)+h*k3');

k1 = u_sum(0,r(i,:));
k2 = u_sum(0,r(i,:)+0.5*h*k1');
k3 = u_sum(0,r(i,:)+0.5*h*k2');
k4 = u_sum(0,r(i,:)+h*k3');

r(i+1,:) = r(i,:) + (h/6)*(k1' + 2*k2' + 2*k3' + k4');

end
toc
plot3(r(:,1),r(:,3),r(:,2),'LineStyle','-','Color','b')
 grid on
 axis([0 2*Lx 0 Lz a b])
 xlabel('x = [0,Lx]'); zlabel('y = [-1,1]'); ylabel('z = [0,Lz]')      % y,z switch to match channelflow convention
 title('Upper Branch Equilibrium')
 hold on

%       end
%   end
 

 plot3(r(size(r,1),1),r(size(r,1),3),r(size(r,1),2),'Marker','.','Color','r')
 plot3(Lx/2,3*Lz/4,0,'LineStyle','-','Marker','.','Color','k')
 hold on
 plot3(Lx,3*Lz/4,0,'LineStyle','-','Marker','.','Color','k')
 