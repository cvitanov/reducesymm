function result = PKxytouvwq(x)


x1 = x(1);
x2 = x(2);
y1 = x(3);
y2 = x(4);

result = [x1^2 + x2^2;
		  y1^2 + y2^2;
		  2*x1^2*y1 + 4*x1*x2*y2 - 2*x2^2*y1;
		  -2*x1^2*y2 + 4*x1*x2*y1 + 2*x2^2*y2];
