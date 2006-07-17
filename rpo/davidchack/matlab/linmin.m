function [p,xi,fret] = linmin(func,p,xi)
% LINMIN - Numerical Recipes routine
%
% [P,XI,FRET] = LINMIN(FUNC,P,XI)
%
%   Given an n-dimensional point P and an n-dimensional direction 
%   XI, moves and resets P to where the function FUNC(P) takes on a 
%   minimum along the direction XI from P, and replaces XI by the 
%   actual vector displacement that P was moved.  Also returns as 
%   FRET the value of FUNC at the returned location P.  This is 
%   actually all accomplished by calling the routines MNBRAK 
%   and BRENT.

% (c) Translated to Matlab by Ruslan L. Davidchack 
%     Date: 20-Nov-2005

  % Global variables communicate with f1dim.
  global pcom xicom nrfunc

  TOL = 2.0e-4; % Tolerance passed to brent.

  % Define the global variables.
  nrfunc = func; pcom = p; xicom = xi;
  ax = 0.0;  xx = 1.0; % Initial guess for brackets.
  [ax,xx,bx,fa,fx,fb] = mnbrak(@f1dim,ax,xx);
%  disp([ax,xx,bx,fa,fx,fb]);
  [fret,xmin] = brent(@f1dim,ax,xx,bx,TOL);
  xi = xmin.*xi;  p = p+xi;
end

function f = f1dim(x)  % Must accompany LINMIN.
  global pcom xicom nrfunc
  f = nrfunc(pcom + x.*xicom);
end
