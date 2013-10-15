function J = Jacobian(t);

%Integrate variational equations before using this function.

load timeevvar.mat;

nt = floor(t/deltat)+1;
n = 4; %number of system dimensions, EDIT THIS

jvec = x(n+1:n^2+n,nt);
J = vector2matrix(jvec);
