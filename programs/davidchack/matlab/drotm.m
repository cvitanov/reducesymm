function [x, jac] = drotm(x0);
% function [x, jac] = drotm(x0) - Double Rotor Map (no wrapping of x(1:2))
%    (optional calculation of Jacobian)

  global L M c

  x = zeros(4,1);
  x(1:2) = M*x0(3:4) + x0(1:2);
  x(3:4) = L*x0(3:4) + c.*sin(x(1:2));

  if nargout > 1,
    jac = [eye(2) M; diag(c.*cos(x(1:2))) zeros(2)];
    jac(3:4,3:4) = L + jac(3:4,1:2)*M;
  end,

return;

%  Example of how to use drotm
np = 20000;
global L M c
f = 8.0; drm_init;
drm = zeros(4,np);
drm(:,1) = [1.0;0.0;0.0;0.1];
for ii = 2:np, drm(:,ii) = drotm(drm(:,ii-1)); end,
plot(drm(3,:),drm(4,:),'r.');
