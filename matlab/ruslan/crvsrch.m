function [sm, xm, hm, nitr] = crvsrch(fcn, xn, crv)
% function [sm, xm, hm, nitr] = crvsrch(fcn, xn, crv)
%
% Given the list of node points XN and the values HN of the scalar function
% HN = H(XN) at these points, CRVSRCH constructs a (rational) polynomial
% interpolation (approximation, extrapolation) curve (selected by input TYPE)
% through XN and finds the point XM which minimizes H(X) along this curve. 
% The function H(X) is the norm of vector function G(X) supplied through FCN.  
% HM = H(XM).

  [nd, np] = size(xn);  % dimension of x and number of node points
% To parametrize the curve c(s), assume points are ordered with 
% c(s=-np+1) = xn(:,1) and c(s=0) = xn(:,np)
  s = (xn-repmat(xn(:,np),1,np))'*(xn(:,1)-xn(:,np))./...
      ((xn(:,1)-xn(:,np))'*(xn(:,1)-xn(:,np))).*(-np+1);
%  disp(s);
  
  global GITR

  switch crv, % types of curve:
  case {1,2,3,4} % polynomial of order crv (np > crv)
    pc = zeros(nd,crv+1);
    for ii = 1:nd, 
%      size(s), size(xn(ii,:)), size(pc(ii,:)),
      pc(ii,:) = polyfit(s',xn(ii,:),crv); end
  case 5, % rational polynomial: order (2,1) (np >= 4)
    nc = zeros(nd,3); dc = zeros(nd,2);
    for ii = 1:nd, [nc(ii,:),dc(ii,:)] = ratpolyfit(s',xn(ii,:),2,1); end
  case 6, % rational polynomial: order (3,1) (np >= 5)
    nc = zeros(nd,4); dc = zeros(nd,2);
    for ii = 1:nd, [nc(ii,:),dc(ii,:)] = ratpolyfit(s',xn(ii,:),3,1); end
  end

  GITR = 0;
  if crv < 5, 
    [as,bs,cs,ha,hb,hc] = mnbrak(@(s)polyf(s,pc,fcn), s(np-1), s(np));
    [hm, sm] = brent(@(s)polyf(s,pc,fcn), as,bs,cs, 1e-2);
    xm = sum(pc'.*repmat(sm.^(size(pc,2)-1:-1:0)',1,size(pc,1)))';
  else, 
    [as,bs,cs,ha,hb,hc] = mnbrak(@(s)ratpolyf(s,nc,dc,fcn), s(np-1), s(np));
    [hm, sm] = brent(@(s)ratpolyf(s,nc,dc,fcn), as,bs,cs, 1e-2); 
    xm = sum(nc'.*repmat(sm.^(size(nc,2)-1:-1:0)',1,size(nc,1)))'./ ...
         sum(dc'.*repmat(sm.^(size(dc,2)-1:-1:0)',1,size(dc,1)))';
  end
  nitr = GITR;
return

function  h = polyf(s, pc, fcn)
  global GITR
  xs = sum(pc'.*repmat(s.^(size(pc,2)-1:-1:0)',1,size(pc,1)))';
  g = fcn(xs);  h = g'*g;  GITR = GITR+1;
return

function h = ratpolyf(s, nc, dc, fcn)
  global GITR
  xs = sum(nc'.*repmat(s.^(size(nc,2)-1:-1:0)',1,size(nc,1)))'./ ...
       sum(dc'.*repmat(s.^(size(dc,2)-1:-1:0)',1,size(dc,1)))';
  g = fcn(xs);  h = g'*g;  GITR = GITR+1;
return