

%open('UB.mat');      % Upper Branch equilibrium
%uu = ans.UB;
open('UB_geom.mat');
ug = ans.UB_geom;

 % geometry settings, don't edit
 Mx = ug(1); My = ug(2); Mz = ug(3); Nx = ug(4); Ny = ug(5); Nz = ug(6);
 Nd = ug(7); Lx = ug(8); Lz = ug(9); a = ug(10); b = ug(11);
 lx = ug(12); lz = ug(13); alpha = ug(14); gamma = ug(15);
 % end geometry settings
 
nx = 143;
ny = 104;
nz = 143;
dx = Lx/nx;
dy = (b-a)/ny;
dz = Lz/nz;

tic
for i = 0:nx
    for j = 0:ny
        for k = 0:nz
            rr = [i*dx a+j*dy k*dz];
            r_grid(1,i+1,j+1,k+1) = rr(1);
            r_grid(2,i+1,j+1,k+1) = rr(2);
            r_grid(3,i+1,j+1,k+1) = rr(3);
        end
    end
end
toc
