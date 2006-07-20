if 0,   %%% Use "global" Newton-Raphson method to improve periodic orbit location
  p = 20;  filename = 'ik16p20a';
  load([filename '.upo']);
  eval(['upos = ' filename ';']); clear(filename);
  xp = upos((1:20)+40,2);  yp = upos((1:20)+40,3);
  
  map = 'ikeda'; par = [1.0;0.9;0.4;6.0]; 
  [xpp,ypp] = feval([map 'xy'],par,p,xp,yp);
  erp = sqrt((xp-xpp).^2+(yp-ypp).^2);

 ip = [4:p 1:3]';
%  xi = xp(ip);  yi = yp(ip);
  xj = mapitr(map, par, [p 0], [xp(ip(1));yp(ip(1))]);
  xi = xj(1,:)';  yi = xj(2,:)';
  for k = 1:20,
    [xx,yy,dxx,dxy,dyx,dyy] = feval([map 'jxy'],par,1,xi,yi);
    g = [[xx(1:p-1)-xi(2:p);xx(p)-xi(1)]...
         [yy(1:p-1)-yi(2:p);yy(p)-yi(1)]]';
    g = g(:);  dg = zeros(2*p);  erg = norm(g);
    for i1 = 1:p-1,
      dg(2*i1+[-1 0],2*i1+[-1:2]) = ...
        [dxx(i1) dxy(i1) -1 0;dyx(i1) dyy(i1) 0 -1];
    end,
    dg(2*p+[-1 0],2*p+[-1 0]) = [dxx(p) dxy(p); dyx(p) dyy(p)];
    dg(2*p+[-1 0],1:2) = -eye(2);
    dx = -inv(dg)*g;  xi = xi + dx(1:2:end);  yi = yi + dx(2:2:end);
    erx = sqrt(sum((xi - xp(ip)).^2 + (yi - yp(ip)).^2));    
    disp(['erg = ' num2str(erg) '   erx = ' num2str(erx)]);
  end,
end,


if 1,   %%% Explore convergence basins for Double Rotor Map
  load dr80ops1-5;  p = 2;  n = 10;
  if ~exist('nsold'), addpath('c:\documents and settings\ruslan l davidchack\work\matlab\newton\solvers\'); end,
  global L M c P
  f = 8.0;  drm_init;  P = p;
  xp = op{p}(n,:)';
  a = (-1:0.01:1)'; na = length(a); randn('seed',123454); r = randn(1,4);  r = r./norm(r); %r = [1 0 0 0 0 0];
  xini = (repmat(xp(:)',na,1) + a*r)';  figure(1);
  for ia = 1:na,
    x = xini(:,ia); g = drotmg(x); 
    plot(a(ia),norm(g),'r*'); hold on;
    if 1, % Newton-Raphson
      x = xini(:,ia);
      for k = 1:10,  [g, jac] = drotmg(x);  dx = -inv(jac)*g;  x = x + dx;  end,
      plot(a(ia),norm(x-xp(:)),'*','color',[0 0 1]); end,
    if 1, % Newton-Armijo 
      x = xini(:,ia);
      [sol, it_hist, ierr, x_hist] = nsold(x,'drotmg',[1e-8; 0]);  x = sol;
      plot(a(ia),norm(x-xp(:))-.1,'*','color',[0 .7 0]);  end,
    if 0, % Newton-Krylov
      x = xini(:,ia);
      [sol, it_hist, ierr, x_hist] = nsoli(x,'ikedaggl',[1e-8; 0]);  x = sol;
      plot(a(ia),norm(x-xp(:))-.1,'*','color',[0 .7 0]);  end,
    if 0, % Broyden
      x = xini(:,ia);
     [sol, it_hist, ierr, x_hist] = brsola(x,'ikedaggl',[1e-8; 0]);  x = sol;
      plot(a(ia),norm(x-xp(:))-.1,'*','color',[0 .7 0]);  end,
    if 1,
      x = xini(:,ia);
      [g, jac] = drotmg(xp(:));  [u,s,C] = svd(jac);  C = -C*u';  beta = 20;
      for k = 1:40,  [g, jac] = drotmg(x);  dx = inv(beta*norm(g)*C' - jac)*g;  x = x + dx;  end,
      plot(a(ia),norm(x-xp(:))-.2,'*','color',[.7 .3 0]);  end,
  end, set(gca,'ylim',[-.5 5]);  
end,


if 0,   %%% Explore convergence basins for global Ikeda map
  load ik16ops15;  p = 6;  n = 2;
  map = 'ikeda';  par = [1.0; 0.9; 0.4; 6.0]; 
  if ~exist('nsold'), addpath('c:\documents and settings\ruslan l davidchack\work\matlab\newton\solvers\'); end,
  global GLOBPAR
  GLOBPAR = par;
  xp = [op(p).x((n-1)*p+(1:p))';  op(p).y((n-1)*p+(1:p))'];
  a = (-3:0.05:3)'; na = length(a); randn('seed',123455); r = randn(1,2*p);  r = r./norm(r); %r = [1 0 0 0 0 0];
  xini = (repmat(xp(:)',na,1) + a*r)';  figure(3);
  for ia = 1:na,
    x = xini(:,ia); g = ikedaggl(x); 
    plot(a(ia),norm(g),'r*'); hold on;
    if 1, % Newton-Raphson
      x = xini(:,ia);
      for k = 1:10,  [g, jac] = ikedaggl(x);  dx = -inv(jac)*g;  x = x + dx;  end,
      plot(a(ia),norm(x-xp(:)),'*','color',[0 0 1]); end,
    if 1, % Newton-Armijo 
      x = xini(:,ia);
      [sol, it_hist, ierr, x_hist] = nsold(x,'ikedaggl',[1e-8; 0]);  x = sol;
      plot(a(ia),norm(x-xp(:))-.1,'*','color',[0 .7 0]);  end,
    if 0, % Newton-Krylov
      x = xini(:,ia);
      [sol, it_hist, ierr, x_hist] = nsoli(x,'ikedaggl',[1e-8; 0]);  x = sol;
      plot(a(ia),norm(x-xp(:))-.1,'*','color',[0 .7 0]);  end,
    if 0, % Broyden
      x = xini(:,ia);
     [sol, it_hist, ierr, x_hist] = brsola(x,'ikedaggl',[1e-8; 0]);  x = sol;
      plot(a(ia),norm(x-xp(:))-.1,'*','color',[0 .7 0]);  end,
    if 1,
      x = xini(:,ia);
      [g, jac] = ikedaggl(xp(:));  [u,s,C] = svd(jac);  C = -C*u';  beta = 10;
      for k = 1:60,  [g, jac] = ikedaggl(x);  dx = inv(beta*norm(g)*C' - jac)*g;  x = x + dx;  end,
      plot(a(ia),norm(x-xp(:))-.2,'*','color',[.7 .3 0]);  end,
  end, set(gca,'ylim',[-.5 5]);  
end,

if 0,    %%% Explore the structure of global map Jacobian and it's stabilisation for UPOs
  load ik16ops15;  p = 3;  n = 1;
  map = 'ikeda';  par = [1.0; 0.9; 0.4; 6.0]; 
  global GLOBPAR;  GLOBPAR = par;
  xp = [op(p).x((n-1)*p+(1:p))';  op(p).y((n-1)*p+(1:p))'];
  [g, dg] = ikedaggl(xp(:));  [c, s, v] = svd(dg);  c = -v*c';
end,

if 0,   %%% Constructing stabilising transformations for arbitrary seeds
  name = 'ikeda';  par = [1.0; 0.9; 0.4; 6.0];  x0 = [0.5; -1.0];  p = 3;
  x = mapitr(name, par, p, x0);  x = x(:);
  [g, dg] = ikedaggl(x);  [u,s,C] = svd(dg);  C = -C*u';   beta = 10;
  for k = 1:200, [g, dg] = ikedaggl(x);  disp(g');  dx = inv(beta*norm(g)*C' - dg)*g;  x = x + dx;  end,
  
end,

if 0,   %%% Exploring Jacobian of Ikeda map
  figure(gcf); clf;
  par = [1.0;0.9;0.4;6.0];
  x = (-.5:0.1:2)'; y = (-2.5:0.1:1)'; [x0,y0] = meshgrid(x,y);
  [xx,yy,dxx,dxy,dyx,dyy] = ikedajxy(par,1,x0,y0);
  subplot(2,3,1); [c,h] = contourf(x,y,xx,[-5 -1 0 1]); clabel(c,h); caxis([-5 5]); axis image;  axis([-.5 2 -2.5 1]); 
  subplot(2,3,4); [c,h] = contourf(x,y,yy,[-5 -1 0 1]); clabel(c,h); caxis([-5 5]);  axis image;  axis([-.5 2 -2.5 1]);
  subplot(2,3,2); [c,h] = contourf(x,y,dxx,[-5 -1 0 1]); clabel(c,h); caxis([-5 5]);  axis image;  axis([-.5 2 -2.5 1]); 
  subplot(2,3,3); [c,h] = contourf(x,y,dxy,[-5 -1 0 1]); clabel(c,h); caxis([-5 5]);  axis image;  axis([-.5 2 -2.5 1]); 
  subplot(2,3,5); [c,h] = contourf(x,y,dyx,[-5 -1 0 1]); clabel(c,h); caxis([-5 5]);  axis image;  axis([-.5 2 -2.5 1]); 
  subplot(2,3,6); [c,h] = contourf(x,y,dyy,[-5 -1 0 1]); clabel(c,h); caxis([-5 5]);  axis image;  axis([-.5 2 -2.5 1]); 
end,