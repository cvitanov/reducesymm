function x1 = ikeda16(x0)
  par = [1.0; 0.9; 0.4; 6.0];
  x = x0(1); y = x0(2);
  [xx,yy] = ikedaxy(par,1,x,y);
  x1 = [xx-x; yy-y];
return
  