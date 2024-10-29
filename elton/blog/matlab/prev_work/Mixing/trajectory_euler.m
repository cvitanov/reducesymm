% trajectory.m
% John Elton, updated April 27, 2008

% Load spectral coefficients for a velocity field and pick an initial
% starting position and this program will evaluate the spectral sum to find
% the velocity field, then repeat and intergrate a trajectory.

u = UB;   % define u to be whatever flowfield was loaded
ug = UB_geom;   %  same for geometry

% geometry settings
Mx = ug(1); My = ug(2); Mz = ug(3); Nx = ug(4); Ny = ug(5); Nz = ug(6);
Nd = ug(7); Lx = ug(8); Lz = ug(9); a = ug(10); b = ug(11);
lx = ug(12); lz = ug(13); alpha = ug(14); gamma = ug(15);

N = size(u,1);

% create a column vector containing the spectral coefficients for each
% component
uspec1 = zeros(N/3,1);
uspec1 = u(1:N/3,1) + i*u(1:N/3,2);

uspec2 = zeros(N/3,1);
uspec2 = u(N/3+1:2*N/3,1) + i*u(N/3+1:2*N/3,2);

uspec3 = zeros(N/3,1);
uspec3 = u(2*N/3+1:N,1) + i*u(2*N/3+1:N,2);

% pick an initial (x,y,z) value
% x(1) = (pi)/(8*1.14); y(1) = cos(5*pi/34); z(1) = (pi)/(10);
% x(1) = (pi)/(3*1.14); y(1) = cos(11*pi/34); z(1) = (2*pi)/(7.5);
x = []; y = []; z = [];
 x(1) = 1/2; y(1) = 1/2; z(1) = 1/2;

%step size and number of integration steps
dt = 0.1;
tf = 200;

% start loop for integration steps
tic
for jj = 1:tf;

zvec = [1/2];
for mz = 1:Mz-1
    zvec = [zvec; exp(2*pi*i*mz*z(jj)/Lz)];
end

xvec = [];
for mx = 0:Mx-1
    if 0 <= mx & mx <= Mx/2
        kx = mx;
    else
        kx = mx - Mx;
    end      
    xvec = [xvec; exp(2*pi*i*kx*x(jj)/Lx)];
end

T = zeros(My,1);                    % should generalize to [a,b].
T(1) = 1;        % 0th cheby poly
T(2) = y(jj);
for my = 3:My
    T(my) = 2*y(jj)*T(my-1) - T(my-2);
end

v = kron(xvec,zvec);
w = kron(T,v);

u1 = 2*real(uspec1.'*w);
u2 = 2*real(uspec2.'*w);
u3 = 2*real(uspec3.'*w);
u = [u1;u2;u3];

% Euler step, for now
x(jj+1) = x(jj) + u1*dt;
y(jj+1) = y(jj) + u2*dt;
z(jj+1) = z(jj) + u3*dt;

end

r = [x' y' z'];

toc


 plot3(x,z,y)
 axis([-Lx Lx -Lz Lz a b])





