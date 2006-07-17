function [ax,bx,cx,fa,fb,fc] = mnbrak(func,ax,bx);
% MNBRAK - Numerical Recipes routine
%
% [AX,BX,CX,FA,FB,FC] = MNBRAK(FUNC,A,B)
%
%   Given a function FUNC, and given distinct initial points A and B, 
%   this routine searches in the downhill direction (defined by the 
%   function as evaluated at the initial points) and returns new 
%   points AX, BX, CX that bracket a minimum of the function.  Also 
%   returned are the function values at the three points, FA, FB, and FC.

% (c) Translated to Matlab by Ruslan L. Davidchack 
%     Date: 19-Nov-2005

  GOLD = 1.618034;  GLIMIT = 100.0;  TINY = 1.0e-20;
% GOLD is the default ratio by which successive intervals are magnified; 
% GLIMIT is the maximum magnification allowed for a parabolic-fit step.

  fa = func(ax);  fb = func(bx);
  if fb > fa, cx = ax;  ax = bx;  bx = cx;  
    fc = fa;  fa = fb;  fb = fc; end
  cx = bx + GOLD.*(bx-ax);  fc = func(cx);
  while fb > fc,
    r = (bx-ax).*(fb-fc);  q = (bx-cx).*(fb-fa);
    u = bx-((bx-cx).*q-(bx-ax).*r)/(2.0*max(abs(q-r),TINY).*sign(q-r));
    ulim = bx+GLIMIT.*(cx-bx);
    if (bx-u).*(u-cx) > 0.0, fu = func(u);
      if fu < fc, ax = bx;  bx = u;  fa = fb;  fb = fu; return;
      elseif fu > fb, cx = u;  fc = fu; return; end
      u = cx+GOLD.*(cx-bx);  fu = func(u);
    elseif (cx-u).*(u-ulim) > 0.0, fu = func(u);
      if fu < fc, bx = cx;  cx = u;  u = cx+GOLD.*(cx-bx);
        fb = fc;  fc = fu;  fu = func(u); end
    elseif (u-ulim).*(ulim-cx) >= 0.0, u = ulim;  fu = func(u);
    else u = cx+GOLD.*(cx-bx);  fu = func(u); end
    ax = bx;  bx = cx;  cx = u;  fa = fb;  fb = fc;  fc = fu;
  end,
end
