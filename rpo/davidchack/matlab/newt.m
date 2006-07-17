function [x,check] = newt(func,x0);
% Given an initial guess x for a root in N dimensions, find the 
% root by a globally convergent Newton’s method. The length N 
% vector of functions to be zeroed, called fvec in the routine
% below, is returned by a user-supplied routine that must be called
% funcv and have the declaration FUNCTION funcv(x). The output 
% quantity CHECK is false on a normal return and true if the 
% routine has converged to a local minimum of the function fmin 
% defined below.  In this case try restarting from a different 
% initial guess.
% Parameters:
%   MAXITS (200) is the maximum number of iterations;
%   TOLF (1e-4) sets the convergence criterion on function values
%   TOLMIN (1e-6) sets the criterion for deciding whether spurious 
%     convergence to a minimum of fmin has occurred; 
%   TOLX (eps) is the convergence criterion on delta_x;
%   STPMX (100.0) is the scaled maximum step length allowed 
%     in line searches.

  global func, fvec;

  MAXITS = 200;  TOLF = 1e-4;  TOLMIN = 1e-6; 
  TOLX = 3e-16;  STPMX=100.0
  
  x = x0(:);  f = fmin(x);  % Also calculates fvec = func(x)
  if (max(abs(fvec)) < 0.01*tolf), % Test for initial guess being
    check = false;                 % a root.  Use more stringent 
    return;  end,                  % test than simply TOLF.

% Calculate stpmax for line searches.
  stpmax = STPMX.*max(norm(x(:)),length(x));
  for its = 1:MAXITS, % Start of iteration loop.
% If analytic Jacobian is available, you can replace the 
% routine fdjac below with your own routine.
    fjac = fdjac(x,fvec);
    g = fvec'*fjac; % Compute grad f = f*Df for the line search.
    xold = x;  fold = f; % Store x and f.
    p = -fjac\fvec; % p = -inv(J)*fvec;
    [x,f,check] = lnsrch(xold,fold,g,p,stpmax,@fmin);
% lnsrch returns new x and f. It also calculates fvec at the 
% new x when it calls fmin.
    if (max(abs(fvec)) < TOLF), % Test for convergence on 
      check = false;            % function values.
      return; end,
    if (check), % Check for gradient of f zero, i.e., spurious convergence.
      check = (max(abs(g).*max(abs(x),1.0)./max(f,0.5*length(x))) < TOLMIN);
      return; end, 
% Test for convergence on dx.
    if (max(abs(x-xold)./max([abs(x);1.0]) < TOLX),
      return; end,
  end do
  disp('MAXITS exceeded in newt');
return;


function [x,f,check] = lnsrch(xold,fold,g,p,stpmax,func);
% function [x,f,check] = lnsrch(xold,fold,g,p,stpmax,func);
% Given an N-dimensional point xold, the value of the function 
% and gradient there, fold and g, and a direction p, finds a new 
% point x along the direction p from xold where the function func 
% has decreased “sufficiently.”  xold, g, p, and x are all vectors 
% of length N.  The new function value is returned in f.  stpmax 
% is an input quantity that limits the length of the steps so that
% you do not try to evaluate the function in regions where it is 
% undefined or subject to overflow.  p is usually the Newton 
% direction. The output quantity check is false on a normal exit. 
% It is true when x is too close to xold or a roundoff problem
% occured.  In a minimization algorithm, this usually signals 
% convergence and can be ignored.  However, in a zero-finding 
% algorithm the calling program should check whether the 
% convergence is spurious.
% Parameters: 
%   ALF (1e-4) ensures sufficient decrease in function value; 
%   TOLX (eps) is the convergence criterion on delta_x.

  ALF = 1e-4;  TOLX = 3e-16;

  check = false;  pabs = norm(p);
% Scale if attempted step is too big.
  if pabs > stpmax, p =p.*(stpmax/pabs); end,
  slope = g'*p;
  if slope >= 0.0, 
    disp('Roundoff problem in lnsrch'); check = true; return; end,
  alamin = TOLX./max(abs(p)./max(abs(xold),1.0)); % Compute lambda_min
  alam = 1.0; % Always try full Newton step first.
  while 1, % Start of iteration loop.
    x = xold+alam*p; f = func(x);
    if alam < alamin, % Convergence on delta_x. For zero finding,
      x = xold;       % the calling program should verify the convergence.
      check = true;  return;
    elseif (f <= fold+ALF*alam*slope), % Sufficient function decrease.
      return;
    else, % Backtrack.
      if alam == 1.0, % First time.
        tmplam = -slope./(2.0*(f-fold-slope));
      else, % Subsequent backtracks.
        rhs1 = f-fold-alam*slope;
        rhs2 = f2-fold-alam2*slope;
        a = (rhs1/alam.^2-rhs2/alam2.^2)/(alam-alam2);
        b = (-alam2*rhs1/alam.^2+alam*rhs2/alam2.^2)/(alam-alam2);
        if a == 0.0, tmplam = -slope/(2.0*b);
        else, disc = b*b-3.0*a*slope;
          if disc < 0.0, tmplam = 0.5*alam;
          elseif b <= 0.0, tmplam = (-b+sqrt(disc))/(3.0*a);
          else tmplam = -slope/(b+sqrt(disc)); end,
        end,
        if tmplam > 0.5*alam, tmplam = 0.5*alam; end, % ? ? 0.5?1.
      end,
    end,
    alam2 = alam;  f2 = f;
    alam = max(tmplam,0.1*alam); % ? ? 0.1?1.
  end, % Try again.
return;


function fjac = fdjac(x, fvec); 
% function fjac = fdjac(x, fvec);

  global func;
  EPS = 1e-7; % EPS is the approximate square root of the machine precision.
  
  n = length(x);
  xsav = x;  h = EPS.*abs(xsav);  h(find(h == 0.0)) = EPS;
  xph = xsav+h;  % Trick to reduce finite precision error.
  h = xph-xsav;
  for j=1:n,
    x(j) = xph(j);
    df(:,j) = (func(x)-fvec)/h(j); % Forward difference formula.
    x(j) = xsav(j);
  end;
return;


function f = fmin(x);
% function f = fmin(x);
  global func, fvec;
  fvec = func(x);  f = 0.5.*norm(fvec);
return;