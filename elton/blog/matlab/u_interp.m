
function u = u_interp(tt,r)

global uu ug Mx My Mz Nx Ny Nz Nd Lx Lz a b lx lz alpha gamma N uspec1 uspec2 uspec3 nx ny nz dx dy dz v_grid r0 jy jz xcord rr

 x = r(1); y = r(2)+1; z = r(3);
 
 if y > 2
     u = [1; 0; 0];
     return
 elseif y < 0
     u = [-1; 0; 0];
     return
 end
% if rr is outside of the periodic box, shift it back inside to calculate the velocity field
if x > Lx 
    div = floor(x/Lx);
    x = x - div*Lx;
elseif x < 0 
    div = floor(abs(x/Lx));
    x = x+(1+div)*Lx;
end
if z > Lz 
    div = floor(z/Lz);
    z = z - div*Lz;
elseif z < 0 
    div = floor(abs(z/Lz));
    z = z+(1+div)*Lz;
end


gx = x/dx;
gy = y/dy;
gz = z/dz;
i = floor(gx); 
j = floor(gy);
k = floor(gz); 
gx = gx-i;
gy = gy-j;
gz = gz-k;
i = i+1;
j = j+1;
k = k+1;
if i>49 || i<1 || k>49 || k<1 || j>36 || j<1
    i
    j
    k
end
if imag(i) ~= 0 || imag(j) ~= 0 || imag(k) ~= 0
    i
    j
    k
    x
    gx
    dx
    r0
    jy 
    jz
    xcord
    rr
end

u1 = (1-gx)*(1-gy)*(1-gz)*v_grid{i,j,k}(1) + gx*(1-gy)*(1-gz)*v_grid{i+1,j,k}(1) + ...
     (1-gx)*gy*(1-gz)*v_grid{i,j+1,k}(1) + gx*gy*(1-gz)*v_grid{i+1,j+1,k}(1) + ...
     (1-gx)*(1-gy)*gz*v_grid{i,j,k+1}(1) + gx*(1-gy)*gz*v_grid{i+1,j,k+1}(1) + ...
     (1-gx)*gy*gz*v_grid{i,j+1,k+1}(1) + gx*gy*gz*v_grid{i+1,j+1,k+1}(1);
     
u2 = (1-gx)*(1-gy)*(1-gz)*v_grid{i,j,k}(2) + gx*(1-gy)*(1-gz)*v_grid{i+1,j,k}(2) + ...
     (1-gx)*gy*(1-gz)*v_grid{i,j+1,k}(2) + gx*gy*(1-gz)*v_grid{i+1,j+1,k}(2) + ...
     (1-gx)*(1-gy)*gz*v_grid{i,j,k+1}(2) + gx*(1-gy)*gz*v_grid{i+1,j,k+1}(2) + ...
     (1-gx)*gy*gz*v_grid{i,j+1,k+1}(2) + gx*gy*gz*v_grid{i+1,j+1,k+1}(2);
     
u3 = (1-gx)*(1-gy)*(1-gz)*v_grid{i,j,k}(3) + gx*(1-gy)*(1-gz)*v_grid{i+1,j,k}(3) + ...
     (1-gx)*gy*(1-gz)*v_grid{i,j+1,k}(3) + gx*gy*(1-gz)*v_grid{i+1,j+1,k}(3) + ...
     (1-gx)*(1-gy)*gz*v_grid{i,j,k+1}(3) + gx*(1-gy)*gz*v_grid{i+1,j,k+1}(3) + ...
     (1-gx)*gy*gz*v_grid{i,j+1,k+1}(3) + gx*gy*gz*v_grid{i+1,j+1,k+1}(3);
     
u = [u1; u2; u3];
     
     

