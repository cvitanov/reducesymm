function [fx,xmin] = brent(f,ax,bx,cx,tol);
% BRENT - Numerical Recipes routine
%
% [FX,Xmin] = BRENT(F,AX,BX,CX,TOL)
%
%   Given a function F and given a bracketing triplet of abscissas 
%   AX, BX, CX [such that BX is between AX and CX, and F(BX) is less 
%   than both F(AX) and F(CX)], this routine isolates the minimum to 
%   a fractional precision of about TOL using Brent’s method.  The 
%   abscissa of the minimum is returned as Xmin, and the minimum 
%   function value is returned as FX.

% (c) Translated to Matlab by Ruslan L. Davidchack 
%     Date: 20-Nov-2005

  ITMAX = 100;  CGOLD = 0.3819660;  ZEPS = 1.0e-18;
  % ITMAX is the maximum allowed number of iterations; 
  % CGOLD is the golden ratio; 
  % ZEPS is a small number that protects against trying to achieve 
  %      fractional accuracy for a minimum that happens to be 
  %      exactly zero.

  a = min(ax,cx);  b = max(ax,cx);
  v = bx;  w = v;  x = v;  e = 0.0;
  fx = f(x);  fv = fx;  fw = fx;
  
%  disp('Iterations in BRENT');
  for iter = 1:ITMAX,
%    disp(sprintf('%3d  %10.5f  %10.5f  %10.2e',iter,x,fx,b-a));
    xm = 0.5*(a+b);  tol1 = tol*abs(x)+ZEPS;  tol2 = 2.0*tol1;
    if abs(x-xm) <= (tol2-0.5*(b-a)),  xmin = x;  return; end
    if abs(e) > tol1,
      r = (x-w).*(fx-fv);  q = (x-v).*(fx-fw);
      p = (x-v).*q-(x-w).*r;  q = 2.0.*(q-r);
      if q > 0.0, p = -p; end,  q = abs(q);
      etemp = e;  e = d;
      if (abs(p) >= abs(0.5.*q.*etemp) | p <= q.*(a-x) | p >= q.*(b-x)),
        if x >= xm, e = a-x; else, e = b-x; end,  d = CGOLD.*e;
      else, d = p/q;  u = x+d; 
        if u-a < tol2 | b-u < tol2, d = tol1.*sign(xm-x); end, end
    else,
      if x >= xm, e = a-x; else, e = b-x; end,  d = CGOLD.*e; end
    if abs(d) >= tol1, u = x+d; else, u = x + tol1.*sign(d); end
    fu = f(u);
    if fu <= fx, if u >= x, a = x; else b = x; end
      v = w;  fv = fw;  w = x;  fw = fx;  x = u;  fx = fu;
    else, if u < x, a = u; else b = u; end
      if fu <= fw | w == x,
        v = w;  fv = fw;  w = u;  fw = fu;
      elseif fu < fv | v == x | v == w,
        v = u;  fv = fu;  end, end, end
  disp('BRENT: Exceed maximum iterations');  fx = NaN;  xmin = NaN;  
end
