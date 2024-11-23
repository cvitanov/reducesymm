
clear
global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N uspec1 uspec2 uspec3 

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
 
nx = 48;
ny = 35;
nz = 48;
dx = Lx/nx;
dy = (b-a)/ny;
dz = Lz/nz;

tic
for i = 0:nx
    for j = 0:ny
        for k = 0:nz
            rr = [i*dx a+j*dy k*dz];
            u = u_sum(0,rr);
            v_g{i+1,j+1,k+1} = u';
        end
    end
end
toc