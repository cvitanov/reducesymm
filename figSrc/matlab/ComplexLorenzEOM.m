function dy = ComplexLorenzEOM(t,y,params)

rho1 = params{1};
rho2 = params{2};
sigma = params{3};
b = params{4};
e = params{5};

x1 = y(1);
x2 = y(2);
y1 = y(3);
y2 = y(4);
z = y(5);

dy = zeros(5,1);

% Complex Lorenz equations
dy(1) = -sigma*x1 + sigma*y1;
dy(2) = -sigma*x2 + sigma*y2;
dy(3) = (rho1 - z)*x1 - rho2*x2 - y1 - e*y2;
dy(4) = rho2*x1 + (rho1 - z)*x2 + e*y1 - y2;
dy(5) = -b*z + x1*y1 + x2*y2;
