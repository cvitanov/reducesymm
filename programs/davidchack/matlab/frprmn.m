function [p,fret,iter] = frprmn(func,p,ftol)
% FRPRMN - Numerical Recipes routine
%
% [P,FRET,ITER] = FRPRMN(FUNC,P,FTOL)
%
%   Given a starting point P, Fletcher-Reeves-Polak-Ribiere 
%   minimization is performed on a function [F,GF] = FUNC(P), 
%   using its gradient GF.  The convergence tolerance on the 
%   function value is input as FTOL.  Returned quantities are
%   P (the location of the minimum), ITER (the number of 
%   iterations that were performed), and FRET (the minimum value 
%   of the function).  The routine LINMIN is called to perform 
%   line minimizations.

  ITMAX = 200;  EPS = 1.0e-18;
  % ITMAX is the maximum allowed number of iterations;
  % EPS is a small number to rectify the special case of 
  % converging to exactly zero function value.

  [fp,xi] = func(p); g = -xi; h = g; xi = h; % Initializations.

  for its = 1:ITMAX, % Loop over iterations.
    iter = its;
    [p,xi,fret] = linmin(func,p,xi);
    % Next statement is the normal return:
    if 2.0*abs(fret-fp) <= ftol.*(abs(fret)+abs(fp)+EPS), 
      return; end,
    fp = fret; [gg,xi] = func(p);  gg = 0.0;  dgg = 0.0;

for (j=1;j<=n;j++) {
    gg += g[j]*g[j];
/* dgg += xi[j]*xi[j]; */ This statement for Fletcher-Reeves.
dgg += (xi[j]+g[j])*xi[j]; This statement for Polak-Ribiere.
}
if (gg == 0.0) { Unlikely. If gradient is exactly zero then
FREEALL we are already done.
return;
}
gam=dgg/gg;
for (j=1;j<=n;j++) {
g[j] = -xi[j];
xi[j]=h[j]=g[j]+gam*h[j];
}
}
nrerror("Too many iterations in frprmn");
}