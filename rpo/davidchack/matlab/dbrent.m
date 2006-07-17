function [fx,xmin] = dbrent(f,ax,bx,cx,tol);
% DBRENT - Numerical Recipes routine
%
% [FX,Xmin] = DBRENT(F,AX,BX,CX,TOL)
%
%   Given a function F and its derivative function DF [both calculated
%   in [FX,DFX] = F(X), where DFX is optional], and given a bracketing 
%   triplet of abscissas AX, BX, CX [such that BX is between AX and CX, 
%   and F(BX) is less than both F(AX) and F(CX)], this routine isolates 
%   the minimum to a fractional precision of about TOL using a 
%   modification of Brent’s method that uses derivatives.  The abscissa 
%   of the minimum is returned as Xmin, and the minimum function value 
%   is returned as FX.

% (c) Translated to Matlab by Ruslan L. Davidchack 
%     Date: 19-Nov-2005

  ITMAX = 100;  ZEPS = 1.0e-18;

  a = min(ax,cx);  b = max(ax,cx);
  v = bx;  w = v;  x = v;  e = 0.0;
  [fx,dx] = f(x);  fv = fx;  fw = fx;  dv = dx;  dw = dx;
  
  disp('Iterations in DBRENT');
  for iter = 1:ITMAX,
    disp(sprintf('%3d  %10.5f  %10.5f  %10.2e',iter,x,fx,b-a));
    xm = 0.5*(a+b);  tol1 = tol*abs(x)+ZEPS;  tol2 = 2.0*tol1;
    if abs(x-xm) <= (tol2-0.5*(b-a)),  xmin = x;  return;  end,
    if abs(e) > tol1, d1 = 2.0*(b-a);  d2 = d1;
      if dw ~= dx, d1 = (w-x)*dx/(dx-dw);  end,
      if dv ~= dx, d2 = (v-x)*dx/(dx-dv);  end,
      u1 = x+d1;  u2 = x+d2;
      ok1 = (a-u1)*(u1-b) > 0.0 & dx*d1 <= 0.0;
      ok2 = (a-u2)*(u2-b) > 0.0 & dx*d2 <= 0.0;
      olde = e;  e = d;
      if ok1 | ok2, 
        if ok1 & ok2, if abs(d1) < abs(d2), d = d1; else d = d2; end,
        elseif ok1, d = d1;
        else d = d2; end,
        if abs(d) <= abs(0.5*olde), u = x+d;
          if u-a < tol2 | b-u < tol2, d = tol1.*sign(xm-x); end,
        else, if dx >= 0.0, e = a-x;  else  e = b-x; end, d = 0.5*e; end,
      else, if dx >= 0.0, e = a-x;  else  e = b-x; end, d = 0.5*e; end,
    else, if dx >= 0.0, e = a-x;  else  e = b-x; end, d = 0.5*e; end,
    if abs(d) >= tol1, u = x+d;  fu = f(u);
    else, u = x + tol1.*sign(d); fu = f(u);
      if fu > fx, xmin = x; return; end, end,
    [fu,du] = f(u);
    if fu <= fx, if u >= x, a = x; else b = x; end,
      v = w;  fv = fw;  dv = dw;
      w = x;  fw = fx;  dw = dx;
      x = u;  fx = fu;  dx = du;
    else, if u < x, a = u; else b = u; end,
      if fu <= fw | w == x,
        v = w;  fv = fw;  dv = dw;
        w = u;  fw = fu;  dw = du;
      elseif fu < fv | v == x | v == w,
        v = u;  fv = fu;  dv = du; end, end, end,
  disp('DBRENT: Exceed maximum iterations');  fx = NaN;  xmin = NaN;  
return;
