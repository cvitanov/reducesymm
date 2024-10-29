keep v_grid 
tic
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

 r0 = [2.35105561774981 0.42293662349708 0.65200166068573];    % sp5
%  r0 = [3.16051044117966 -0.42293662349708 0.60463540075018];   % sp6
%  r0 = [0 0 3*Lz/4]    % sp1: [Lx/2 0 Lz/4]  sp2: [Lx/2 0 3*Lz/4] sp3: [0 0 Lz/4]  sp4: [0 0 3*Lz/4]
 
 u_sum(0,r0);
 
%  tic
%  ui = u_interp(0,r0)
% toc


