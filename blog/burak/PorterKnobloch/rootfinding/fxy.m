function result = fxy(x, u)

%Roots of this function is the x1,x2,y1,y2 corresponding to invariant 
%polynomials u,v,w,q

x1 = x(1);
x2 = x(2);
y1 = x(3);
y2 = x(4);

u = u(1);
v = u(2);
w = u(3);
q = u(4);

result = [x1^2 + x2^2 - u;
		  y1^2 + y2^2 - v;
		  2*x1^2*y1 + 4*x1*x2*y2 - 2*x2^2*y1 - w;
		  -2*x1^2*y2 + 4*x1*x2*y1 + 2*x2^2*y2 - q];
