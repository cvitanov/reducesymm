function dudt = kszolflow(u, dx, sig);
  c0 = 0.5./dx;  c1 = dx.^(-2).*(2.0-6.0.*dx.^(-2));
  c2 = dx.^(-2).*(4.0.*dx.^(-2)-1.0);  c3 = dx.^(-4);

  dudt = zeros(size(u));  n = length(u);
  dudt(1) = (c1-c0.*u(2)).*u(1) + c2.*u(2) - c3.*(u(3)+sig.*u(1));
  dudt(2) = (c1-c0.*(u(3)-u(1))).*u(2) + c2.*(u(3)+u(1)) - c3.*u(4);
  dudt(3:n-2) = (c1-c0.*(u(4:n-1)-u(2:n-3))).*u(3:n-2) + ...
    c2.*(u(4:n-1)+u(2:n-3)) - c3.*(u(5:n)+u(1:n-4));
  dudt(n-1) = (c1-c0.*(u(n)-u(n-2))).*u(n-1) + c2.*(u(n)+u(n-2)) - ...
    c3.*u(n-3);
  dudt(n) = (c1+c0.*u(n-1)).*u(n) + c2.*u(n-1) - ...
    c3.*(sig.*u(n)+u(n-2));
