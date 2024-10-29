
%global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N

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
 
 % create a 1-d array containing the spectral coefficients, for each component, with real/im parts interlaced  
uut = uu'; % now the real/im parts are alternating in memory (in matlab, the inner most dimension changes fastest in memory )
uspc1 = uut(:,1:N/3); % real/im parts interlaced in memory for uspec1 
uspc2 = uut(:,N/3+1:2*N/3);
uspc3 = uut(:,2*N/3+1:N);
 
 % end geometry settings
 
 ti = 100;               % starting time
 tf = 0;              % final time
 h = -.1;                % step size
 np = (tf-ti)/h;      % number of steps
 rmax = 0.001; 
 
 % compute some paramters used in c runge program
 % Put the stuff needed for the c runge program that does not change into one array
 % except keep uspecs separate
 parms = [h np  Lx Lz Mx My Mz];
 
%  for jy = 1:Py-1
%      xcord = rmax - 2*rmax*(jy-1)/(Py-2);
%      rr = sqrt(rmax^2 - xcord^2);
%      for jz = 1:Pz-1        
%         r0 = [Lx/2+xcord rr*sin(2*pi*(jz-1)/(Py-2)) Lz/4+rr*cos(2*pi*(jz-1)/(Py-2))];
 
         r0 = -(1/10)*evecs_eq2(:,1);
                  r0 = [r0(1)+Lx/2 r0(2) r0(3)+3*Lz/4]; 
tic
r = runge_sum_c(r0,parms,uspc1,uspc2,uspc3); % will return a 3x(np+1) array 
toc

plot3(r(1,:),r(3,:),r(2,:),'LineStyle','-','Color','b')
 grid on
 axis([0 2*Lx 0 Lz a b])
 xlabel('x = [0,Lx]'); zlabel('y = [-1,1]'); ylabel('z = [0,Lz]')      % y,z switch to match channelflow convention
 title('Upper Branch Equilibrium')
 hold on

%       end
%   end
 

 plot3(r(1,size(r,2)),r(3,size(r,2)),r(2,size(r,2)),'Marker','.','Color','r')
 plot3(Lx/2,3*Lz/4,0,'LineStyle','-','Marker','.','Color','k')
 hold on
 plot3(Lx,3*Lz/4,0,'LineStyle','-','Marker','.','Color','k')
 