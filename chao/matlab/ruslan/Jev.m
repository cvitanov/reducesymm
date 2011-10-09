function [ ev ] = Jev( J, ev0, idx )
% Floquet eigenvectors transport along ppo.
% Given Floquet eigenvectors ev0 of a pre-periodic orbit at a point x_0,
% compute them at points x(t_idx) along the periodic orbit by
% multiplying with 'partial' Jacobian J^{t_idx}(x_0). Here J stores 
% multiple partial Jacobians along the orbit.
% ES, 2011-10-09.

N=size(J,1);
JJ=J(:,(idx-1)*N+1:idx*N);

ev=JJ*ev0;

end

