keep v_grid sp5 sp6
global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N uspec1 uspec2 uspec3 dx dy dz  

open('UB.mat');      % Upper Branch equilibrium
uu = ans.UB;
open('UB_geom.mat');
ug = ans.UB_geom;

 % geometry settings, don't edit
 Mx = ug(1); My = ug(2); Mz = ug(3); Nx = ug(4); Ny = ug(5); Nz = ug(6);
 Nd = ug(7); Lx = ug(8); Lz = ug(9); a = ug(10); b = ug(11);
 lx = ug(12); lz = ug(13); alpha = ug(14); gamma = ug(15);
 N = size(uu,1);
nx = 144;
ny = 105;
nz = 144;
dx = Lx/nx;
dy = (b-a)/ny;
dz = Lz/nz;

eps = 5*10^-4;
u2small = 1;

tic
for i = 0:nx-1
    for j = 0:ny-1
        for k = 0:nz-1
           u = v_grid{i+1,j+1,k+1};
           u2 = u(1).^2 + u(2).^2 + u(3).^2;
           r = [i*dx a+j*dy k*dz];
           if u2 < eps % & 2<r(1) & 3>r(1) & 0.2<r(2) & 0.6>r(2) & 0.5<r(3) & 1>r(3)
               plot3(r(1),r(3),r(2),'Marker','.','Color','b')
               grid on
               axis([0 2*Lx 0 Lz a b])
               xlabel('x = [0,2Lx]'); zlabel('y = [-1,1]'); ylabel('z = [0,Lz]')      % y,z switch to match channelflow convention
               title('Upper Branch Equilibrium')
               hold on 
%                if u2 < u2small
%                    u2small = u2;
%                    rsmall = r;
%                end
           end
        end
    end
end

%plot3(rsmall,'Marker','*','Color','k')
 plot3(Lx/2,3*Lz/4,0,'LineStyle','-','Marker','.','Color','k')
 plot3(Lx/2,Lz/4,0,'LineStyle','-','Marker','.','Color','k')
 plot3(Lx,3*Lz/4,0,'LineStyle','-','Marker','.','Color','k')
 plot3(Lx,Lz/4,0,'LineStyle','-','Marker','.','Color','k')
 plot3(0,3*Lz/4,0,'LineStyle','-','Marker','.','Color','k')
 plot3(0,Lz/4,0,'LineStyle','-','Marker','.','Color','k')
 plot3(sp5(1),sp5(3),sp5(2),'LineStyle','-','Marker','.','Color','k')
 plot3(sp6(1),sp6(3),sp6(2),'LineStyle','-','Marker','.','Color','k')
 
 
toc