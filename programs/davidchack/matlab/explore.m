if 1,   %%% Check DFT of u^2
  uk = [0:7] + [0 1i*(7:-1:1)]; uk = [uk 0 fliplr(conj(uk(2:end)))]';
  u = real(ifft(uk));  uk2 = fft(u.^2).*16;
  iu = mod(repmat((0:15)',1,16)-repmat(0:15,16,1),16)+1;
  uk22 = uk(iu)*uk;
end,


if 0,   %%% Test basins on the Ikeda map
  clear;  par = [1.0; 0.9; 0.4; 6.0];  p = 3;
  load ik16ops15;  xp = [op(1).x; op(p).x];  yp = [op(1).y; op(p).y];  np = length(xp);
  bndr = [-1.0 2.0 -2.8 1.1];
  [x0, y0] = bndrgrid(bndr,400000);
  x1 = x0(:);  y1 = y0(:);  ipc = zeros(length(x1),2);  i0 = (1:length(x1))';  disp(length(i0));
  if 1,
    c11 = zeros(length(i0),1);  c12 = c11;  c21 = c11;  c22 = c11;
    [gx,gy] = iikedaxy(par, p, x1(i0), y1(i0));
    [gx,gy,gxx,gxy,gyx,gyy] = ikedajxy(par,p,gx,gy);
    for ii = i0',
      [u,s] = schur([gxx(ii) gxy(ii); gyx(ii) gyy(ii)]);
%      [u,s,v] = svd([gxx(ii) gxy(ii); gyx(ii) gyy(ii)]);
      c = eye(2) - 2.*u(:,1)*u(:,1)';
%      [u,s,v] = svd(u*s*u'-eye(2));  c = -v*u';
      c11(ii) = c(1,1);  c12(ii) = c(1,2);  c21(ii) = c(2,1);   c22(ii) = c(2,2);
    end,
  end,
  
  for itr = 1:80,
    [gx,gy,gxx,gxy,gyx,gyy] = ikedajxy(par,p,x1(i0),y1(i0));
    gx = gx - x1(i0);  gy = gy - y1(i0);  gxx = gxx - 1;  gyy = gyy - 1;
  
    switch 3,
    case 1,   %  Newton-Raphson
      ng = gxx.*gyy-gyx.*gxy;  dx =  (gxy.*gy-gyy.*gx)./ng;  dy = (-gxx.*gy+gyx.*gx)./ng;
    case 2,   %  Newton-Armijo
      ng = gxx.*gyy-gyx.*gxy;  dx =  (gxy.*gy-gyy.*gx)./ng;  dy = (-gxx.*gy+gyx.*gx)./ng;
      ng = sqrt(gx.^2+gy.^2);  ii = (1:length(dx))';
      for jj = 0:8,
        xa = x1(i0(ii))+dx(ii)./(2^jj);  ya = y1(i0(ii))+dy(ii)./(2^jj);  
        [gx,gy] = ikedaxy(par, p, xa, ya);
        iin = find(sqrt((gx-xa).^2+(gy-ya).^2) < ng(ii).*(1-2.^(-jj-1)));
        if ~isempty(iin), dx(ii(iin)) = dx(ii(iin))./(2^jj);
          dy(ii(iin)) = dy(ii(iin))./(2^jj);  ii(iin) = [];  end,
      end,
      if ~isempty(ii), dx(ii) = 3.*bndr(2) - 2.*bndr(1);  ii = [];  end,
    case 3,   %  DL with stabilising matrix C
    % c11 = 1;  c12 = 0;  c21 = 0;  c22 = 1;
    % c11 = -0.76225788674941;  c12 = 0.64727344614801;  c21 = 0.64727344614801;  c22 = 0.76225788674941;
    % c11 = -0.48253549731930;  c12 = 0.87587641469948;  c21 = 0.87587641469948;  c22 = 0.48253549731930;
    % c11 = -0.14158995133880;  c12 =-0.98992539399688;  c21 =-0.98992539399688;  c22 = 0.14158995133880;
    % c11 = 0;  c12 = -1;  c21 = -1;  c22 = 0;
      
      beta = 20;  ng = sqrt(gx.^2+gy.^2);  b = beta.*ng;
      ng = b.^2.*(c11(i0).*c22(i0)-c12(i0).*c21(i0)) + gxx.*gyy-gyx.*gxy + b.*(c12(i0).*gxy+c21(i0).*gyx-c11(i0).*gyy-c22(i0).*gxx);
      dx = (b.*(c22(i0).*gx-c21(i0).*gy)+gxy.*gy-gyy.*gx)./ng;
      dy = (b.*(c11(i0).*gy-c12(i0).*gx)+gyx.*gx-gxx.*gy)./ng;
    end,
    
    x1(i0) = x1(i0) + dx;  y1(i0) = y1(i0) + dy;
    
    ii = find(x1(i0) < bndr(1) | x1(i0) > bndr(2) | y1(i0) < bndr(3) | y1(i0) > bndr(4));
    if ~isempty(ii),  ipc(i0(ii),1) = 0;  ipc(i0(ii),2) = itr;  iout = ii;  end,
    for ip = 1:np,
      ii = find(abs(x1(i0)-xp(ip))<1e-8 & abs(y1(i0)-yp(ip))<1e-8);
      if ~isempty(ii), ipc(i0(ii),1) = ip; ipc(i0(ii),2) = itr;  iout = [iout; ii];  
        disp(sprintf('%3d  %3d  %10d',itr,ip,length(ii)));  end,
    end,
    if ~isempty(iout), i0(iout) = [];  iout = [];  end, disp(sprintf('%3d  %10d',itr,length(i0)));  
    if isempty(i0), return; end,
  end,
  ipc(i0,1) = 0;  ipc(i0,2) = itr;
  
  if 1,  figure(1); clf; cmap = hsv(8); cmap = [[1 1 1];cmap]; 
    h = plotcolr(x0(:),y0(:),ipc(:,1)+1,cmap); set(h,'markersize',2); axis image;
    hold on;  plot(xp,yp,'k*');  hold off;  set(gcf,'pos',[0 60 720 885]);
  end,

  return;
end,


if 0,   %%% Test methods for integrating KSFM flow
  N = 32;  nu = 0.015;  k = (1:N)';  J = k.^2.*(1-nu.*k.^2);
%  randn('seed',1230000); a0 = randn(N,1).*exp(-0.3*k);  
  load ksfm7seed;   
  b = a;  h = 2e-4;  figure(1); hold on;  hpl = plot(a(1),a(2),'r-','markersize',2); grid on;
%  Wi = 1./(1./h - J);  %%% case 2
  gam = 1./(2+sqrt(2));  Wi = 1./(1 - h.*gam.*J);  %%% case 3
  for ii = 1:14400,
    switch 4, 
    case 1,  %%% Euler
      a = a + h.*ksfmflow(a);
    case 2,  %%% W-method (order 1) 
      a = a + Wi.*ksfmflow(a);
    case 3,  %%% W-method (order 2)
      k1 = Wi.*ksfmflow(a);
      k2 = Wi.*(ksfmflow(a+0.5.*h.*k1)-k1)+k1;
      a = a + h.*k2;  %if mod(ii,2)==0, b = [b a];  end,
    case 4,
      k1 = Wi.*ksfmflow(a);
      k2 = Wi.*(ksfmflow(a+(2.0./3.0).*h.*k1)-(4.0/3.0).*h.*gam.*J.*k1);
      a = a + h.*(k1 + 3.0.*k2)./4.0;
    end,
    set(hpl,'xdata',[get(hpl,'xdata') a(1)],'ydata',[get(hpl,'ydata') a(2)]);
    if mod(ii,2000)==0 | abs(a(1)) > 10, disp(h*ii); pause, end,
  end,
  return;
end,


if 0,  %%% Test analytic Jacobian calculation
  eps = 1e-6;  load ksfmseed11;  global C;  C = eye(N);
  y = [a; ti];
  if 1,  dg = ksfmc3j(ti, y);  
  else   [g, dg] = ksfmflowmapj(nu, ti, a, 2e-4);  end,
  ddg = zeros(size(dg));
  for ii = 1:N+1,
    if 1, 
      yp = y;  yp(ii) = yp(ii) + eps;  ym = y;  ym(ii) = ym(ii) - eps;
      gp = ksfmc3(ti, yp);  gm = ksfmc3(ti, ym);
    else,  ap = a;  ap(ii) = ap(ii) + eps;  am = a;  am(ii) = am(ii) - eps;
      gp = ksfmflowmapj(nu, ti, ap, 2e-4);  gm = ksfmflowmapj(nu, ti, am, 2e-4);  end,
    ddg(:, ii) = (gp-gm)./(2.*eps);
  end,
  figure(1); clf; mesh(dg-ddg); return;
end,
if 0,
  eps = 1e-5;  load ksfmseed11;  
  [f, df] = ksfmj(nu, a);  ddf = zeros(size(df));
  for ii = 1:N,
    ap = a;  ap(ii) = ap(ii) + eps;  am = a;  am(ii) = am(ii) - eps;
    fp = ksfmj(nu, ap);  fm = ksfmj(nu, am);
    ddf(:,ii) = (fp-fm)./(2.*eps);
  end,
  figure(1); clf; mesh(df-ddf); return;
end,


if 0,  %%% SD iteration scheme with explicit time variable
  if 0,
    N = 32;  nu = 0.015;
    ai = zeros(N,1);  ai(1) = -1.27;  ti = 1.0;
    tic; [b,t] = flowdp8('ksfm', [N; nu], ai, 0, 0:0.02:ti, 3e-4); toc; 
    x = (0:.01:1)'*pi; u = zeros(length(x),length(t));
    for i=1:length(t), u(:,i) = sum(repmat(b(:,i),1,length(x)).*sin(x*(1:N))')'; end,
    figure(1); clf; pcolor(x,t,u'); shading flat; pause;
    a = b(:,end);  f = ksfm(nu,a);  %disp(norm(f)); f = f./norm(f);
    ti = 3.0;
    [b,s] = flowdp8('ksfm', [N; nu], a, 0, 0:0.02:ti, 3e-4);
    t = 1.44;  np = 5;
    figure(1); clf;
    global C h1 h2 h3; C = eye(N);
    h1 = plot(s,sum(repmat(f./norm(f),1,size(b,2)).*(b-repmat(a,1,size(b,2)))),'b.-'); hold on; grid on;
    h2 = plot(s,sqrt(sum((repmat(a,1,size(b,2))-b).^2)),'.-','color',[0 .8 0]);
    figure(2); clf;
    h3 = plot(b(1,:),b(2,:),'.-'); hold on; grid on;
  else  load ksfm7seed;  end,
  if 1,
    N = 32;  nu = 0.015;  h = 2e-4;
    if 1,
      apre = zeros(N,1);  apre(1) = -2.23;  
      randn('seed',12340011); apre = randn(N,1).*exp(-0.3.*(1:N)');
      tpre = 1.3;  ti = 3.0;
      [a, df] = ksfmflowmapj(nu, tpre, apre, h);  f = ksfm(nu, a);
    else,
      load ksfmseed12;  f = ksfm(nu, a);
    end,
    s = (0:0.01:1)'*ti;
    b = ksfmflowmapj(nu, s, a, h);
    global C h1 h2 h3;
    [u,r,v] = svd((eye(N)-0.*f*f')*df);  nud = 3;  nst = 1;  disp([find(diag(r)>1)' nst]);  disp(duposp(nud,nst)); 
    C = eye(N) + u(:,1:nud)*(duposp(nud,nst)-eye(nud))*u(:,1:nud)';
    if 1,
      dd = b-repmat(a,1,size(b,2));
      figure(1);  clf;  h1 = plot(s,sum(repmat(f./norm(f),1,size(b,2)).*dd),'b.-'); hold on; grid on;
      h2 = plot(s,sqrt(sum(dd.^2)),'.-','color',[0 .8 0]);
      figure(2); clf;  h3 = plot(b(1,:),b(2,:),'.-'); hold on; grid on;
    end,
    t = 1.0;%  return;
  end,
  if 0,  
    load ksfm147upo;  
    global C;  [y,r] = schur(jac);  C = eye(N) - 2.*(y(:,1)*y(:,1)');
    a = a0; a(1) = a(1) + 1e-6;  t = tend;
  end,
  if 1,  %%% Use Matlab solvers
    reltol = 1e-6;  abstol = reltol;
    options = odeset('abstol',abstol,'reltol',reltol); %,'jacobian',@ksfmc4j);
    [s,y] = ode15s(@ksfmc4, [0 100], [a; t],options);
  return; end,
  if 0,  %%% Use radau5Mex solver
    opt.RelTol = 1e-5*ones(33,1);  opt.AbsTol = opt.RelTol(1).*[10*exp(-0.3*(1:32)'); 1];
    opt.Jacobian = 'ksfmc3j';  % opt.JacobianLowerBandwidth = 0;  
    % opt.JacobianUpperBandwidth = 0;  opt.InitialStep = 0.005;
    [s,y] = radau5Mex('ksfmc4', [0 10.0], [a; t], opt);
  return; end,
  h = 1e-2;  tt = 0;  nacc = 0;  so = 0;  soo = 0;  nreg = 0;  clear do go dd;
  ii = 0;  nitr = input('NITR = ');
  while 1,
    ii = ii + 1;
    switch 7,
    case 1,  %%% Euler steps without PSS
      [b,s] = flowdp8('ksfm', [N; nu], a, 0, [0:0.02:1]*t, 2e-4);
      ft = ksfm(nu,b(:,end));  ft = ft./norm(ft);
      g = b(:,end)-a;
      if exist('go'), dd = g - go;  dd = dd./norm(dd); else err = 0; end,
      if exist('do'), err = 1 - do'*dd; else err = 0; end,
      disp(sprintf('%3d  %9.5f  %9.5f  %9.5f  %9.5f  %9.5f  %9.5f  %9.5f  %2d  %2d',ii,norm(g),t,err,h,tt,so,soo,nacc,nreg)); pause;
      if err < 0.5,  % accept the step
        pln = sum(repmat(ft,1,size(b,2)).*(b-repmat(a,1,size(b,2))));
        pld = sqrt(sum((repmat(a,1,size(b,2))-b).^2));
        figure(1);  set(h1,'xdata',s,'ydata',pln);  plot(s(end),pln(end),'ko');
        set(h2,'xdata',s,'ydata',pld); plot(s(end),pld(end),'ro');
        figure(2);  set(h3,'xdata',b(1,:),'ydata',b(2,:));
        plot(b(1,1),b(2,1),'ro',b(1,end),b(2,end),'ko');
        if nacc > 3 & err < 0.02, h = 1.5*h; end,
        if exist('do'), soo = so;  aoo = ao; too = to; goo = go; foo = fo; doo = do;  end,
        so = tt;  ao = a;  to = t;  go = g;  fo = ft;  if exist('dd'), do = dd; end,
        a = a + h*g;  t = t - h.*(ft'*g);  tt = tt + h;
        nacc = nacc + 1;  nreg = 0;
      else
        nacc = 0;
        h = 0.6.*h;  nreg = nreg + 1;
        if nreg < 2,
          a = ao + h.*go;  t = to - h.*(fo'*go);  tt = so + h;
        else,
        %  so = soo;  ao = aoo;  to = too;  go = goo;  fo = foo;  do = doo;
        %  a = ao + h.*go;  t = to - h.*(fo'*go);  tt = so + h;
        tt = .5*(so+soo);  t = .5*(to+too);  a = .5*(ao + aoo);  clear go do dd;
        end,
      end,
    case 2,  %%% Euler steps with PSS
      [b,s] = pssitr('ksfm', [N; nu], np, a, 0, 2e-4);
      if ~exist('sold'), sold = s(end); end,
      if (s(end) - sold < -0.03),
        np = np + 1;  [b,s] = pssitr('ksfm', [N; nu], np, a, 0, 2e-4);
      elseif (s(end) - sold > 0.03),
        np = np - 1;  b(:,end) = [];   s(end) = [];
      end,
      if (abs(s(end)-sold) > 0.01), cnt = input('Continue? '); end,
      b = [a b];  s = [0.0 s];  sold = s(end);
      ft = ksfm(nu,b(:,end));  ft = ft./norm(ft);
      dn = sum(repmat(ft,1,size(b,2)).*(b-repmat(a,1,size(b,2))));
      dd = sqrt(sum((repmat(a,1,size(b,2))-b).^2));
      figure(1);  set(h1,'xdata',s,'ydata',dn);  plot(s(end),dn(end),'ko');
      set(h2,'xdata',s,'ydata',dd); plot(s(end),dd(end),'ro');
      figure(2);
      set(h3,'xdata',b(1,:),'ydata',b(2,:));
      plot(b(1,1),b(2,1),'ro',b(1,end),b(2,end),'ko');
      g = b(:,end)-a;  disp([ii norm(g) s(end)]);  drawnow;
      g = g - 2.*y(:,1)*(y(:,1)'*g);  %%%
      fg = f'*g;
      a = a + 0.0005*(g - fg.*f);
      f = ksfm(nu,a);  f = f./norm(f);
    case 3,  %%% DL iterations
      b0 = [a eye(N)];  beta = 1e3;
      [b,s] = flowdp8('ksfmt', [N*(N+1); nu], b0(:), 0, [0:0.1:1]*t, 2e-4);
      g = b(1:N,end)-a;  ng = norm(g);  df = reshape(b(N+1:end,end),N,N);
      f = ksfm(nu,a);  ft = ksfm(nu,b(1:N,end));
      dg = [df-eye(N) ft; f' 0]; ft = ft./norm(ft);

      dn = sum(repmat(ft,1,size(b,2)).*(b(1:N,:)-repmat(a,1,size(b,2))));
      dd = sqrt(sum((repmat(a,1,size(b,2))-b(1:N,:)).^2));
      figure(1);  set(h1,'xdata',s,'ydata',dn);  plot(s(end),dn(end),'ko');
      set(h2,'xdata',s,'ydata',dd); plot(s(end),dd(end),'ro');
      figure(2);  set(h3,'xdata',b(1,:),'ydata',b(2,:));
      plot(b(1,1),b(2,1),'ro',b(1,end),b(2,end),'ko');
      disp([ii ng t]); drawnow;
      
      at = [a; t] + (beta*ng*eye(N+1) - dg)\[g; -(ft'*g)];
      a = at(1:N); t = at(N+1);
    case 4,  %%% Adaptive stepsize control with RK1-2 methods
      y0 = [a; t];
      if ii > 1, if err < 1,  g0 = g12;  end,
      else g0 = ksfmgyj2(y0);  hd = 0; hdd = 0; end,
      y11 = y0 + h*g0;             g11 = ksfmgyj2(y11);
      y12 = y0 + 0.5*h.*(g0+g11);  g12 = ksfmgyj2(y12);
      
      erry = 2.*norm(y12-y11)./(norm(y12-y0)+norm(y11-y0));
      errg = 2.*norm(g12-g11)./(norm(g12-g0)+norm(y11-y0));
%      plot(log10(h), log10([erry errg]), '*'); hold on;  h = h*sqrt(2); grid on; err = 2; pause; 
      
      toly = 0.2;  tolg = 0.7;  safety = 0.8;
      err = max([erry./toly errg./tolg]);
      if err < 1.0,
        a = y12(1:N);  t = y12(N+1);
        disp(sprintf('%5d %11.7f %9.2e %9.2e %9.2f  %9.3f', ii, h, erry, errg, errg./erry, abs(hdd/hd)));
      else, disp(sprintf('%5d %11.7f %9.2e %9.2e %9.2f  %9.3f - reject', ii, h, erry, errg, errg./erry, abs(hdd/hd)));  end,
      h = safety/err*h;  hd = 0.2*(safety/err-1) + 0.8*hd;  hdd = 0.2*abs(safety/err-1) + 0.8*hdd;
      figure(3); plot(y0(1),y0(2),'ro', g0(1)+y0(1),g0(2)+y0(2),'bo'); hold on;  grid on;  drawnow;
    case 5,  %%% Adaptive stepsize control with Ros1-2 methods
      y0 = [a; t];  gam = sqrt(0.5);
      [g0, jac] = ksfmgyj(y0);  
      y11 = y0 + (eye(N+1)./h - jac)\g0;  g11 = ksfmgyj(y11);
      k1 = (eye(N+1)./h - gam.*jac)\g0;  g1 = ksfmgyj(y0 + k1);
      y12 = y0 + 1.5.*k1 + (eye(N+1)./h - gam.*jac)\(0.5.*g1 - k1./h);
      g12 = ksfmgyj(y12);
       
      erry = 2.*norm(y11-y12)./(norm(y11-y0)+norm(y12-y0));
      errg = 2.*norm(g11-g12)./(norm(g11-g0)+norm(g12-g0));
      
      toly = 0.2;  tolg = 0.7;  safety = 0.8;
      if erry < toly & errg < tolg,  
        a = y12(1:N);  t = y12(N+1);
        disp(sprintf('%5d %11.7f %9.5f %9.5f %9.2f', ii, h, erry./toly, errg./tolg, errg./erry));
      else, disp(sprintf('%5d %11.7f %9.5f %9.5f %9.2f - reject', ii, h, erry./toly, errg./tolg, errg./erry));  end,
      h = safety*min([toly./erry tolg./errg])*h;
      figure(3); plot(y0(1),y0(2),'m*', g0(1)+y0(1),g0(2)+y0(2),'c*'); hold on;  grid on;  drawnow;
    case 6,   %%% Adaptive 2-step control
      y0 = [a; t];  g0 = ksfmgyj2(y0);
      y1 = y0 + h.*g0;  g1 = ksfmgyj2(y1);
      y2 = y1 + h.*g1;  g2 = ksfmgyj2(y2);
      
      erry = norm(y2-2.*y1+y0)./norm(y2-y0);
      errg = norm(g2-2.*g1+g0)./norm(g2-g0);
      
      toly = 0.1;  tolg = 0.3;  safety = 0.9;
      err = max([erry./toly errg./tolg]);   if ii == 1, erro = err;  end,
      hnew = h.*min([5.0 erro.^-.2.*err.^-.8]).*safety;
      if err < 1,  
        a = y2(1:N);  t = y2(N+1);
        disp(sprintf('%5d %11.7f %9.5f %9.5f %9.2f', ii, h, erry./toly, errg./tolg, errg./erry));
      else, disp(sprintf('%5d %11.7f %9.5f %9.5f %9.2f - reject', ii, h, erry./toly, errg./tolg, errg./erry));  end,
      h = hnew;  erro = err;
      figure(3); plot(y0(1),y0(2),'ro', g0(1)+y0(1),g0(2)+y0(2),'bo'); hold on;  grid on;  pause;
    case 7,   %%% Adaptive RK5 scheme (Dopri5 or HIHA5)
      y0 = [a; t];  g0 = ksfmgyj2(y0);
      [ynew, gnew, erry, errg] = dop5step(@ksfmgyj2, y0, g0, h);
%      disp(sprintf('%5d %11.7f %9.2e %9.2e %9.2f', ii, h, erry, errg, errg./erry));

      toly = 1e-6;  tolg = 2e-6;  err = max([erry./toly errg./tolg]);  safety = 0.8;  
      if err < 1,  
        a = ynew(1:N);  t = ynew(N+1);
        disp(sprintf('%5d %11.7f %9.2e %9.2e %9.2f', ii, h, erry, errg, errg./erry));
        plot(ynew(1),ynew(2),'ro', gnew(1)+ynew(1),gnew(2)+ynew(2),'bo');  hold on;  grid on;  drawnow;
      else, disp(sprintf('%5d %11.7f %9.2e %9.2e %9.2f - reject', ii, h, erry, errg, errg./erry));
        plot(ynew(1),ynew(2),'rx', gnew(1)+ynew(1),gnew(2)+ynew(2),'bx');  hold on;  grid on;  drawnow;
      end,
      h = safety.*h/err^0.25;
 %     plot(log10(h), log10([erry errg]), '*'); hold on;  h = h*sqrt(2);  pause;
    end,
    if ii == nitr,  nitr = nitr + input('NITR = ');  end,
  end,
  return;
end,


if 0,   %%% Test dopri5 and radau5 mex drivers with ksfm
    N = 32;  nu = 0.015;
    ai = zeros(N,1);  ai(1) = -1.27;  ti = 10.0;
    opt.RelTol = 1e-9;  opt.AbsTol = 1e-10;
    opt.Jacobian = 'ksfmj_mex';  opt.JacobianLowerBandwidth = 0;  
    opt.JacobianUpperBandwidth = 0;  opt.InitialStep = 0.0;
    tic; [t,a,stats,tnext] = radau5Mex('ksfm_mex', [0 ti], ai, opt); toc;
    x = (0:.01:1)'*pi; u = zeros(length(x),length(t));
    for i=1:length(t), u(:,i) = sum(repmat(a(i,:)',1,length(x)).*sin(x*(1:N))')'; end,
    figure(2); clf; pcolor(x,t,u'); shading flat;
end,

if 0,   %%% Locating periodic orbits in KS (Fourier modes) with PS adot(1) = 0
  name = 'ksfm';  N = 32;  nu = 0.015;  par = [N; nu];
  a0 = zeros(par(1),1);  a0(1) = -1;  t0 = 0;  h = 2e-4;
  [a,t] = pssitr(name, par, 1, a0, t0, h); disp(t); a0 = a;
  plot(a0(1),a0(2),'*'); hold on;
  [a,t] = pssitr(name, par, 1, a0, t0, h); disp(t);
  plot(a(1),a(2),'r*');
  g = a - a0;
  [f, pss, dps] = ksfm(nu, a);  P = f*dps'./(f'*dps);

  name = 'ksfmt';  par(1) = N*(N+1);
  b0 = [a0 eye(N)];  b0 = b0(:);
  [b,t] = pssitr(name, par, 1, b0, t0, h); disp(t);
  Df = (eye(N) - P)*reshape(b(N+1:end),N,N);
  
  a1 = a0 - 0.01*((Df-eye(N))\g);
  plot(a1(1),a1(2),'g*');
%  [f, pss, dps] = ksfm(nu, a1);  disp(pss);
  
  if 0,T
    for i = 1:4,
      [f, pss, dps] = ksfm(nu, a1);  disp(pss);
      a1 = a1 - pss*dps./norm(dps).^2;  end,
  else
    name = 'ksfm';  par(1) = N;
    [a,t] = pssitr(name, par, 1, a1, t0, h); disp(t);
    a1 = a;  end,
  plot(a1(1),a1(2),'c*');  grid on;
  name = 'ksfm';  par(1) = N;
  [a,t] = pssitr(name, par, 2, a1, t0, h); disp(t);
  plot(a(1,2),a(2,2),'m*');  
end,

if 0,   %%% Locating periodic orbits in KS (Fourier modes) - Locating close approaches
  name = 'ksfm';  N = 32;  nu = 0.015;  par = [N; nu];
  a0 = zeros(par(1),1);  a0(1) = -1.0;  t0 = 0;  tspan = 1.02:0.02:7; h = 2e-4;
  [a,t] = flowdp8(name, par, a0, t0, tspan, h);
 
  x = (0:.01:1)'*pi; u = zeros(length(x),length(t));
  for i=1:length(t), u(:,i) = sum(repmat(a(:,i),1,length(x)).*sin(x*(1:par(1)))')'; end,
  figure(1); clf; pcolor(x,t,u'); shading flat;  %pause;

  if 0, cr = zeros(size(a,2));
    for ii = 1:size(a,2),  cr(ii,:) = sqrt(sum((a-repmat(a(:,ii),1,size(a,2))).^2));  end,
    figure(2); clf; pcolor(cr); shading flat;  pause;
  end,
 
  isd = [48 81];  tend = 0.65914234;
%  isd = [83 103]; tend = 0.396;
%  isd = [224 249];

%  isd = [62 122];
%  isd = [208 216]; tend = 0.161;
  figure(1); hold on; plot(repmat([0 pi]',1,2), repmat(t(isd),2,1),'k-');
  
  if 1, % Explore converfence of svd iterates
    name = 'ksfmt';  par = [N*(N+1); nu];  h = 2e-4;  figure(2); clf;  U = [];
    for ii = 1:40,
      b0 = [a(:,48-ii) eye(N)];  b0 = b0(:);  tend = 0.02*ii;
      [b,t] = flowdp8(name, par, b0, 0.0, [0 .5 1]*tend, h);
      df1 = reshape(b(N+1:end,3),N,N); [u1,s1,v] = svd(df1);
      U(:,:,ii) = u1;
      clf; pcolor(u1); axis image; title(num2str(ii)); caxis([-1 1]); colorbar; pause;
    end,
    return;
  end,

  if 0, % Construct stabilising transformations based on svd of iterates
   name = 'ksfmt';  par = [N*(N+1); nu];
   b0 = [a(:,48) eye(N)];  b0 = b0(:);  tend = 0.65914234;
%   b0 = [a(:,14) eye(N)];  b0 = b0(:);  tend = 0.68;
   [b,t] = flowdp8(name, par, b0, 0.0, [0 .5 1]*tend, h);
   df1 = reshape(b(N+1:end,3),N,N); [u1,s1,v] = svd(df1);
   q1 = u1*(s1-eye(N))*u1';
   s1(1,1) = -s1(1,1);   q2=u1*(s1-eye(N))*u1';
   s1(2,2) = -s1(2,2);   q3=u1*(s1-eye(N))*u1';
   s1(1,1) = -s1(1,1);   q4=u1*(s1-eye(N))*u1';
   [u,s,v] = svd(q1); c1 = -v*u';  [u,s,v] = svd(q2); c2 = -v*u';
   [u,s,v] = svd(q3); c3 = -v*u';  [u,s,v] = svd(q4); c4 = -v*u';
%   b0 = [a(:,2) eye(N)];  b0 = b0(:);  tend = 0.92;
%   b0 = [a(:,48) eye(N)];  b0 = b0(:);  tend = 0.65914234;
%   [b,t] = flowdp8(name, par, b0, 0.0, [0 .5 1]*tend, h);
%   df2 = reshape(b(N+1:end,3),N,N); [u2,s2,v] = svd(df2);
%   q1 = u1*(s1-eye(N))*u1';  q2 = u2*(s2-eye(N))*u2';
%   [u,s,v] = svd(q1);  c1 = -v*u';  [u,s,v] = svd(q2);  c2 = -v*u';
   return;
  end
  
  if 0, 
    a0 = a(:,isd(1));%  tend = 0.02*(diff(isd)+2);
%    [a,t] = flowdp8(name, par, a0, t0, [0 tend], h);
    figure(3); clf; plot(t,sqrt(sum((a-repmat(a0,1,size(a,2))).^2)),'.-'); grid on;
    f0 = ksfm(nu,a0); ps = f0'*(a-repmat(a0,1,size(a,2)))./norm(f0); hold on;  plot(t,ps,'r.-'); return;
  end,
  
end,


if 0, % Explore difference between eigenvectors, svd, schur for pre- and post-iterates
%   clear;  load ksfm2seed;
   b0 = [a0 eye(N)];  b0 = b0(:);  tend = 0.3871;
   [b,t] = flowdp8('ksfmt', [N*(N+1); nu], b0, 0.0, [0 .5 1]*tend, 2e-4);
   df0 = reshape(b(N+1:end,3),N,N);  disp(norm(b(1:N,3)-a0));
%   load ksfmupo2; 
%   b1 = [a00 eye(N)];  b1 = b1(:);%  tend = 0.65914234;
%   [b,t] = flowdp8('ksfmt', par, b1, 0.0, [0 .5 1]*tend, h);
%   df1 = reshape(b(N+1:end,3),N,N);  disp(norm(b(1:N,3)-a00));
   [p0,e0] = eig(df0);     %  [p1,e1] = eig(df1);
   [y0,t0] = schur(df0);   %  [y1,t1] = schur(df1);
   [u0,s0,v0] = svd(df0);  %  [u1,s1,v1] = svd(df1);
   return;
end

  
if 0,  % SD iterations with subspace projection (on PSS)
%  clear;  load ksfmU;
%  a00 = a(:,48);  
  load ksfm1seed; ksi = y0(:,1:2);  C = [0 1; 1 0];
  name = 'ksfm';  par = [N; nu];  h = 2e-4;%fi  a0 = a00; % figure(1); %clf;
  [a,t] = pssitr(name, par, 2, a0, 0, h);  g0 = a(:,2)-a0;  g0 = g0./norm(g0);
  f0 = ksfm(nu,a0);
  for ii = 1:100,
    [a,t] = pssitr(name, par, 2, a0, 0, h);
    g = a(:,2)-a0; % disp([g0'*g f0'*g]./norm(g));
    b = ksi'*g;
    disp([norm(g) t(2)]);
    plot(a0(1),a0(2),'*','color',[1 0 0]); hold on; pause;
    a0 = a0 + 0.005*(ksi*C*b + g - ksi*b);
  end,
  return;
end,


if 0,  % SD iterations (without PSS - explicit time)
  clear;  load ksfmU;
  a00 = a(:,48);  t00 = 0.65914234;  %ksi = U(:,1:2,40);  C = [0 1; -1 0];
  name = 'ksfm';  par = [N; nu];  h = 2e-4;  a0 = a00;  t0 = t00;  figure(1); %clf;
  f00 = ksfm(nu,a00);
  for ii = 1:100,
    f0 = ksfm(nu,a0);
    [a,t] = flowdp8(name, par, a0, 0, [0 .5 1]*t0, h);
    ga = (a(:,3)-a0); % g = g./norm(g); disp([g0'*g f0'*g]);
    gt = f0'*(ga);
    disp([a0(1) norm(a(:,3)-a0) t0 gt]);
    plot(a0(1),a0(2),'*','color',[.6 0 .6]); hold on; pause;
    a0 = a0 + 0.005*ga;  t0 = t0 + 0.005*gt;
  end,
  return;
end,


if 0,  %%% Locating periodic orbits in KS (Fourier modes) - Zoldi's approach
  
%  clear;  load ksfm8seed;  a0 = a;  tend = t;  N = 32; nu = 0.015;
%  load ksfmU;  a0 = a(:,83);  tend = 0.38;%  f0 = ksfm(nu,a0);
  N = 32; nu = 0.015;  a0 = y(end,1:N)';  tend = y(end,N+1);
%  load ksfmseed11;  a0 = a;  tend = ti;
  for itr = 1:10,
    name = 'ksfmt';  par = [N*(N+1); nu];  h = 2e-4;
%    b0 = [a0 eye(N)];  b0 = b0(:);
%    [b,t] = flowdp8(name, par, b0, 0.0, [0 .5 1]*tend, h);
  
%    g = a0 - b(1:N,3);  ng = norm(g);  disp([tend ng]);
%    Dg = reshape(b(N+1:end,3),N,N) - eye(N);
  [b, jac] = ksfmflowmapj(nu, tend, a0, h);  g = a0 - b;  Dg = jac - eye(N);  ng = norm(g);  
  disp(sprintf('%8.5f  %11.3e', tend, ng));
  if 0,
    name = 'ksfm';  par = [N; nu];  dd = zeros(N,1);  dd(1) = 1e-6;
    [ap,t] = flowdp8(name, par, a0+dd, t0, [0 .5 1]*tend, h);
    [am,t] = flowdp8(name, par, a0-dd, t0, [0 .5 1]*tend, h);
    df1 = (ap(:,3)-am(:,3))./2e-6;
    return;
  end,
  
    f0 = ksfm(nu,a0);  
    ft = ksfm(nu,b(1:N));
    invDg = inv([Dg ft; f0' 0]);
    
    for ii = 0:10,  % Armijo iteration
      da = 2.^(-ii).*(invDg*[g; 0]);
      a1 = a0+da(1:N);  te1 = tend + da(N+1);
      if te1 < 2.0*tend & te1 > 0.5*tend,
        name = 'ksfm';  par = [N; nu];
%        [a,t] = flowdp8(name, par, a1, 0.0, [0 .5 1]*te1, h);
        a = ksfmflowmapj(nu, te1, a1, h);
        g1 = a1 - a;  disp(sprintf('    %2d  %11.3e  %11.3e',ii,norm(da),norm(g1)));
        if norm(g1) < ng.*(1-2.^(-ii-2)),  a0 = a1;  tend = te1;  break;  end
        
      end,
    end, pause;
  end,
end,


if 0,   %%% Test if the stationary point of finite difference KS equation is real
  clear;  load ks49odd;  L = (par(1)+1)*par(2);
  if 1,
    for ii = 1:4,
      [f, df] = ksfd(0, x0, L, -1);  disp(norm(f));
      x0 = x0 - (df\f);  end,
  end,
  plot((1:length(x0)).*L./(length(x0)+1), [x0 f], 'b.-'); hold on;
%  t0 = 0;  tspan = 0:50;
%  [x,t] = flowdp8(name, par, x0, t0, tspan, 0.01);
%  figure(1); clf; pcolor((1:par(1)).*par(2),t,x'); shading flat;
  y0 = interp([0; x0; 0],2);  y0 = y0(2:100);
  for ii = 1:6,
    [f, df] = ksfd(0, y0, L, -1);  disp(norm(f));
    y0 = y0 - (df\f);  end,
  plot((1:length(y0)).*L./(length(y0)+1), [y0 f], 'r.-'); hold on;  
  z0 = interp([0; y0; 0],2);  z0 = z0(2:200);
  for ii = 1:6,
    [f, df] = ksfd(0, z0, L, -1);  disp(norm(f));
    z0 = z0 - (df\f);  end,
  plot((1:length(z0)).*L./(length(z0)+1), [z0 f], 'k.-'); hold on;  
%  par = [99; 0.25];
%  [y,t] = flowdp8(name, par, y0, t0, tspan, 0.002);
%  figure(2); clf; pcolor((1:par(1)).*par(2),t,y'); shading flat;
end,

if 0,   %%% Check symmetry of periodic orbits in CHM3 (Three coupled Henon maps)
  clear;
  name = 'cph3';  p = 10;  eps = 1e-6;
  nameupo = [name sprintf('p%02d',p)];
  load([nameupo '.upo']);
  eval(['ip = find(' nameupo '(:,1) == ' num2str(p) ');']);  np = length(ip);
  eval(['upos = ' nameupo '(ip,[2 4 6]);']);
  nper = zeros(np,6);
  for ip = 1:np,
    nper(ip,1) = ip;
    for pp = 2:6,
      per = perm(3,pp)*[1 2 3]';
      ipp = find(abs(upos(ip,per(1))-upos(:,1))<eps & abs(upos(ip,per(2))-upos(:,2))<eps & ...
            abs(upos(ip,per(3))-upos(:,3))<eps);
      if length(ipp) == 1,
        nper(ip,pp) = ipp;
      elseif length(ipp) == 0,
        disp('Missing orbit point:'); disp([ip pp]); disp([upos(ip,per)]);  nper(ip,pp) = 0;
      else
        disp(['ERROR: more than one orbit point match']); disp([ip pp]); disp(ipp);  return;  end,  
    end,
  end,
end,    
    
if 0,   %%% Check symmetry of periodic orbits in DRM
  load dr80ops6;  p = 6;  np = size(op{p},1);  upos = op{p};
  for ip = 1:np,
     rupo = [mod(-upos(ip,1:2)+0.5*pi,2*pi)-0.5*pi -upos(ip,3:4)];
     ir = find(abs(upos(:,1)-rupo(1)) < 1e-8 & abs(upos(:,2)-rupo(2)) < 1e-8 & ...
          abs(upos(:,3)-rupo(3)) < 1e-8 & abs(upos(:,4)-rupo(4)) < 1e-8);
     if (length(ir) ~= 1),  disp([p ip]);  end,
  end,
end,


if 0,   %%% Determine unstable direction for inverse iterates (Ikeda map)
  par = [1.0;0.9;0.4;6.0]; scale = 0.01;
  load ik16ops15;
%  x0 = op(1).x(1);  y0 = op(1).y(1);
  x0 = 0;  y0 = -0.0;
  for ip = 1:2,
    [x y] = iikedaxy(par,ip,x0,y0);
    xp = x + scale*cos((0:.05:2).*pi);  yp = y + scale*sin((0:.05:2).*pi);
    [x0 y0 F11 F12 F21 F22] = ikedajxy(par,ip,x,y);  F = [F11 F12; F21 F22];
    [P,E] = eig(F);
    [U,S,V] = svd(F);
    [xx yy] = ikedaxy(par,ip,xp,yp);
    dy = F*[xp-x; yp-y];
    plot(x0,y0,'k*',xx,yy,'b-',x0+dy(1,:),y0+dy(2,:),'r-'); hold on;
    pe = 0.01*S*U;
    quiver(x0,y0,scale*S(1,1)*U(1,1),scale*S(1,1)*U(1,2));
  end,
  axis image;
end,


if 0,   %%% Investigate construction of stabilising transformations without periodic orbits in DRM
  load dr80ops1-5;  p = 3; 
  
  np = size(op{p},1);
  for ip = 1:np,
    df = jp{p}(:,:,ip)';  e = eig(df);
    [v,s,u] = svd(df);  ns = length(find(s > 1));
    q1 = v*(s-eye(4))*v';  
    s(1,1) = -s(1,1);   q2=v*(s-eye(4))*v';
    s(2,2) = -s(2,2);   q3=v*(s-eye(4))*v';
    s(1,1) = -s(1,1);   q4=v*(s-eye(4))*v';
    [u,s,v] = svd(q1); c1 = -v*u';  [u,s,v] = svd(q2); c2 = -v*u';
    [u,s,v] = svd(q3); c3 = -v*u';  [u,s,v] = svd(q4); c4 = -v*u';
    [e eig(c1*(df-eye(4))) eig(c2*(df-eye(4))) eig(c3*(df-eye(4))) eig(c4*(df-eye(4)))], pause;
%    if min(eig(df)) > -1, col = 'r-';  else col = 'b-'; end,
%    [u,s,v] = svd(df - eye(2));  c = -v*u';
%    plot(x0(ip)+[c(1,1) 0 c(1,2)]*.05,y0(ip)+[c(2,1) 0 c(2,2)]*.05,col);  hold on;
 end,  
% axis equal;
end,


if 0,   %%% Investigate construction of stabilising transformations without periodic orbits in Ikeda map
  load ik16ops15;  p = 7;  np = length(op(p).x);
  par = [1.0;0.9;0.4;6.0];  x0 = op(p).x;  y0 = op(p).y;
  [xx,yy,dxx,dxy,dyx,dyy] = ikedajxy(par,p,x0,y0);
  clf; 
  for ip = 1:np,
    df = [dxx(ip) dxy(ip); dyx(ip) dyy(ip)];
    [u,s,v] = svd(df);  q1 = v*(s-eye(2))*v';  
    s(1,1) = -s(1,1);   q2=v*(s-eye(2))*v';
    [u,s,v] = svd(q1); c1 = -v*u';  [u,s,v] = svd(q2); c2 = -v*u';
    [eig(df) eig(c1*(df-eye(2))) eig(c2*(df-eye(2)))], pause
%    if min(eig(df)) > -1, col = 'r-';  else col = 'b-'; end,
%    [u,s,v] = svd(df - eye(2));  c = -v*u';
%    plot(x0(ip)+[c(1,1) 0 c(1,2)]*.05,y0(ip)+[c(2,1) 0 c(2,2)]*.05,col);  hold on;
 end,
 axis equal;
end,


if 0,   % Iterates of DKP on a grid
  name = 'dkp';  par = -0.1;  bndr = [-.2 .2 -.2 .2];  ytx = 1;
  [n,pn] = bndrgrid(bndr,200000,ytx);  ng = length(n(:));  n1 = n;  pn1 = pn;
 
  for ig = 1:ng,  
    n0 = n(ig);  pn0 = pn(ig);  m0 = 0;
    pm0 = sqrt(4 + 2*par*n0^2 - pn0^2);
    if imag(pm0) == 0,
      x0 = [n0; pn0; m0; pm0];  t0 = 0;  tol = 1e-10;
      [x,t] = pssitr(name, par, [1 1], x0, t0, 0.0, tol);
      n1(ig) = x(1,1); pn1(ig) = x(2,1);
    else
      n1(ig) = NaN; pn1(ig) = NaN; 
    end,
  end,
% contour(n,pn,n1-n,[0 0],'b'); hold on; contour(n,pn,pn1-pn,[0 0],'r');
  figure(1);  clf;
  h = plotcol(n,pn,n1-n,cmap,[-.1 .1]); set(h,'markersize',2);pixgrid(size(n)); axis(bndr);
  figure(2); clf;
  h = plotcol(n,pn,pn1-pn,cmap,[-.1 .1]); set(h,'markersize',2);ixgrid(size(n)); axis(bndr);
end,


if 0,   % Lyapunov spectrum of the coupled Henon maps
  name = 'hencp';  a = 1.4;  b = 0.3;  L = 4;  
  epss = 0:0.0002:0.25;
  mlexp = [];  slexp = [];
  figure(1); clf;
  for eps = epss,
    par = [a; b; L; eps];  x0 = 0.2*rand(2*L,1);
    lexp = lexpmap(name, par, [10000 40000 10], x0);
    errorbar(eps*ones(1,2*L),mean(lexp),std(lexp),'o'); set(gca,'xlim',[0 .25]); hold on; drawnow;
    mlexp = [mlexp mean(lexp)'];  slexp = [slexp 2*std(lexp)'];  disp([eps mean(lexp)]);
  end,
  save hencp4 a b L epss mlexp slexp
end,


if 0,   % Lyapunov spectrum for KS (finite differences)
  load ks01lexp;
  name = 'ks01';  N = 49;  dxs = 0.405:0.01:0.6;
  % mlexp = [];  slexp = [];
  for dx = dxs,
    nu = (pi/dx/(N+1)).^2;
    par = [N; dx; -1];
    x0 = randn(par(1),1); h = 0.01;
    lexp = lexpflow(name, par, [2000 1000 7 50], x0, h, 5);
    errorbar(nu*ones(1,5),mean(lexp),std(lexp),'o'); set(gca,'xlim',[0.01 .025]); hold on; drawnow;
    mlexp = [mlexp mean(lexp)'];  slexp = [slexp 2*std(lexp)'];  disp([nu mean(lexp)]);
  end,
  save ks01lexp name N dxs mlexp slexp
end,


if 0,   % Lyapunov spectrum for KS (Fourier modes)
 name = 'ksfm';  d = 32;  lpar = [];  mlexp = [];  slexp = [];
% load ksfmlexp;
 for nu = 0.011:0.001:0.029,
   par = [d; nu];
   a0 = zeros(par(1),1);  a0(1) = -1;  a0(2:5) = 0.1*randn(4,1);  h = 2e-4;
   lexp = lexpflow(name, par, [2000 1000 7 50], a0, h, 5);
   errorbar(nu*ones(1,5),mean(lexp),std(lexp),'o'); set(gca,'xlim',[0.01 .025]); hold on; drawnow;
   lpar = [lpar nu];  mlexp = [mlexp mean(lexp)'];  slexp = [slexp 2*std(lexp)'];  disp([nu mean(lexp)]);
  end,
  save ksfmlexp name d lpar mlexp slexp;
end,


if 0,   % Perturbation of symmetric matrix by orthogonal matrix 
%  clear;  load dr80ops1-5;  p = 4;  ip = 100;
%  Df = jp{p}(:,:,ip)';
  n = 5;  figure(1);  clf;
  Df = randn(n); G = Df - eye(n);  [u,s,v] = svd(G);  C = v*u';
  for im = 1:100,
    nnp = 0;
    for nn = 1:10000,  P = eye(n);
      for ii = 1:n-1,
        for jj = ii+1:n,
          ang = .7*pi*(rand(1)-0.5);  Q = eye(n);  
          Q(ii,ii) = cos(ang);   Q(ii,jj) = sin(ang);
          Q(jj,ii) = -Q(ii,jj);  Q(jj,jj) = Q(ii,ii); P = Q*P;
        end,  end,
      miep = min(real(eig(P)));  mieg = min(real(eig(P*C*G)));
      if miep > 0, nnp = nnp + 1; end,  plot(miep,mieg,'.','markersize',2); hold on;
      if miep > 0 & mieg < 0, disp(nn); disp(C*G); disp(P);  return; end,
    end, grid on;  disp(nnp); drawnow;
  end,
end,


if 0,   % Connection between eigensystem and polar decomposition
  clear;  syms a b c d s w al th ux vx uy vy
  F = [a b;c d];  G = F - eye(2);  C = [s*cos(al) sin(al);-s*sin(al) cos(al)];
%  E = eig(C*G);
  E = [s w;-w s];
%  E = [s 0;0 w];
  V = [cos(a) exp(b); sin(a) 0];
%  V = [ux vx; uy vy];
  G = V*(E-eye(2))*inv(V);
  [rs1,how] = simple(G(1,2)-G(2,1));
  [rc1,how] = simple(-G(1,1)-G(2,2));
  [rs2,how] = simple(-G(1,2)-G(2,1));
  [rc2,how] = simple(G(1,1)-G(2,2));
end,


if 0,   % UPOs eigenvector directions of the double rotor map
  clear;  load dr80ops1-5;  p = 4;
  ip = find(abs(op{p}(:,3)) < 21.5 & abs(op{p}(:,4)) < 21.5);  np = length(ip);
  figure(1); clf;
  plot(op{p}(ip,1),op{p}(ip,2),'o','markersize',4); hold on;
  vv = [];
  for ii = ip',
    df = jp{p}(:,:,ii)';  
    [v,e] = eig(df);%  disp(ii);  disp(diag(e)');
%      if (sum(abs(imag(diag(e)))) ~= 0),    disp(sprintf('%4d ',[ip(ii) diag(e)'])); end,
      [m,im] = max(abs(diag(e))); if imag(e(im,im)) ~= 0, disp(ii);  disp(diag(e)'); end,
    vv = [vv; v(1:2,im)'];
  end,
  plot([1; 1]*op{p}(ip,1)'+[-1; 1]*vv(:,1)',[1;1]*op{p}(ip,2)'+[-1; 1]*vv(:,2)','r-');
end,


if 0,   % Influence of perturbation by orthogonal matrix on eigenvalues
  clear;  load dr80ops1-5;  p = 4;  ip = 112;
  Df = jp{p}(:,:,ip)';  G = Df - eye(4);  [u,s,v] = svd(G);  C = -v*u';
  figure(1); clf;
  for nn = 1:10,
    A = randn(4);  [u,s,v] = svd(A);
    for th = (-1:0.01:1).*pi,
      P = v*[cos(th) sin(th) 0 0;-sin(th) cos(th) 0 0;0 0 1 0; 0 0 0 1]*v'; ei = eig(P*C*G);
      plot(th/pi*2,real(ei),'.');  hold on;
    end,
  end,
  grid on;
end,


if 0,  % Stabilising with polar Q of seed orbits for double rotor map
  clear;  load dr80ops1-5;  ps = 3;  p = 4;  nps = size(jp{ps},3);   np = size(jp{p},3); 
  for ips = 1:1,
    Df = jp{ps}(:,:,ips)';  [V,E] = eig(Df);  iV = inv(V);
    for ic = 1:4,
      G = V*E*duposp(4,ic-1)*iV - eye(4);  [u,s,v] = svd(G);  C = -v*u';
      ns{ic} = 0;  is{ic} = [];
      for ip = 1:np,
        se = eig(C*(jp{p}(:,:,ip)'-eye(4)));
        if max(se) < 0, ns{ic} = ns{ic} + 1;  is{ic} = [is{ic}; ip];  end,
      end,
    end,
  end,
end,


if 0,  % Transforming orthogonal matrix into SD matrix
  j = [-1 -1 -4; 4  3 -3; -4 -5 -3];
  [u,s,v] = svd(j);
  q = -v*u';
  c = -[0 0 -1; 1 0 0; 0 -1 0]';
  p = c*q';
end,


if 0,  % Calculating rotation numbers for UPOs of the DRM
  clear;  load dr80ops7;
  
  global L M c;  f = 8.0;  drm_init;
  for ip = 7:7,
    np = size(op{ip},1);
    rn{ip} = zeros(np,2);
    for ii = 1:np,
      x = op{ip}(ii,:)';
      for jj = 1:ip,
        x = drotm(x);
      end,
      rn{ip}(ii,:) = round((x(1:2)' - op{ip}(ii,1:2))./(2.*pi));
%      disp(max(abs(x' - [2*pi*rn{ip}(ii,:) 0 0] - op{ip}(ii,:))));
    end,
    nrn{ip} = [];  trn = rn{ip};
    while ~isempty(trn),
      rn1 = trn(1,:);
      irn = find(trn(:,1) == rn1(1) & trn(:,2) == rn1(2));
      nrn{ip} = [nrn{ip}; [rn1 length(irn)./ip]];
      trn(irn,:) = [];
    end,
  end,
  save dr80ops7 op rn nrn
end,


if 0,  % Plot numbers of UPOs with given rotation number
  figure(1); clf;
  clear;  load dr80ops7;  ip = 7;
  plot(repmat([-8.5; 8.5],1,18),repmat(-8.5:8.5,2,1),'k:',repmat(-8.5:8.5,2,1),repmat([-8.5; 8.5],1,18),'k:'); 
  text(nrn{ip}(:,1),nrn{ip}(:,2)-nrn{ip}(:,1),num2str(nrn{ip}(:,3)),'hor','c','ver','m');
  axis equal;  set(gca,'xlim',[-8.5 8.5],'ylim',[-8.5 8.5],'xtick',-8:8,'ytick',-8:8);
  title(['DRM:  Number of UPOs with rotation numbers (N_1,N_2):  Period ' num2str(ip)],'fontname','times','fontsize',12);
  xlabel('N_1','fontname','times','fontsize',12);  ylabel('N_2 - N_1','fontname','times','fontsize',12);
  orient tall;
end,
