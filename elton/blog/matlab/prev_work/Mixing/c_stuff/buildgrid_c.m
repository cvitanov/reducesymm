
%global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N uspec1 uspec2 uspec3 v_grid r_grid dx dy dz

open('UB.mat');      % Upper Branch equilibrium
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
 
nx = 143;
ny = 104;
nz = 143;
% v_grid and r_grid will have dimensions (3,nx+1,ny+1,nz+1)
dx = Lx/nx;
dy = (b-a)/ny;
dz = Lz/nz;

tic
[v_grid,r_grid] = mk_vgrid_c([nx ny nz dx dy dz Mx My Mz Lx Lz],uspc1,uspc2,uspc3);
%[v_grid1,r_grid1] = mk_vgrid_c([nx ny nz dx1 dy1 dz1 Mx My Mz Lx Lz],uspec1_re,uspec1_im,uspec2_re,uspec2_im,uspec3_re,uspec3_im);
% for i = 0:nx
%     for j = 0:ny
%         for k = 0:nz
%             rr = [i*dx a+j*dy k*dz];
%             u = u_sum(0,rr);
%             v_grid(1,i+1,j+1,k+1) = u(1);
%             v_grid(2,i+1,j+1,k+1) = u(2);
%             v_grid(3,i+1,j+1,k+1) = u(3);
%             r_grid(1,i+1,j+1,k+1) = rr(1);
%             r_grid(2,i+1,j+1,k+1) = rr(2);
%             r_grid(3,i+1,j+1,k+1) = rr(3);
%         end
%     end
% end
toc
