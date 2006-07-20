if 1,  %%% Locate close returns for Zoldi's KS
  N = 45;  dx = 0.5;  sig = -1;  h = 0.02; randn('seed',12340000);
  x = (1:N).*dx; u0 = 0.2*randn(N,1);  tspan = 0:500;
  [u,t] = flowdp8('ks01', [N dx sig], u0, 0.0, tspan, h);
  figure(1); clf;  set(gcf,'pos',[5  300  670  640]);
  pcolor(x,t,u'); shading flat;
  
  nt = size(u,2);  cr = zeros(nt);
  for ii = 1:nt,
    cr(ii,:) = sqrt(sum((u-repmat(u(:,ii),1,nt)).^2));
    cr(ii,:) = [cr(ii,ii:end) zeros(1,ii-1)];  end,
  if 1,  figure(2); clf; set(gcf,'pos',[690 300 700 650]); 
    contourf(t,t,cr',1:3:16);  grid on; 
    xlabel('t_0');  ylabel('T','rotat',0);  axis equal;
    [nb, ne] = ginput(1);  ne = ne + nb;
    nb = round(nb./(t(2)-t(1)));  ne = round(ne./(t(2)-t(1)));
    disp(sprintf('    nb = %4d;  ne = %4d; %%',nb,ne));
  else
%  N = 99;  dx = 0.5;  sig = 1;  randn('seed',12340000);
    nb =   89;  ne =  140; %
%    nb =  142;  ne =  168; %
%  N = 99;  dx = 0.5;  sig = 1;  randn('seed',12340001);
%    nb =  150;  ne =  169; %
%  N = 99;  dx = 0.5;  sig = 1;  h = 0.02; randn('seed',12340002);
%    nb =   82;  ne =  112; % T = 36.03: st - 2600, na - 21000
%  N = 99;  dx = 0.5;  sig = 1;  h = 0.02; randn('seed',12340004);
%    nb =  147;  ne =  157; % T = 0 (steady state): st - 700 
%  N = 99;  dx = 0.5;  sig = 1;  h = 0.02; randn('seed',12340012);
    nb =   98;  ne =  165; % close return;
  end,
  figure(1); hold on; plot(x([1 end end 1]),t([nb nb ne ne]),'k-');
  figure(2); clf;  plot(x,u(:,nb),'ro-',x,u(:,ne),'ko-');  grid on;
  figure(3); clf;  
  hds = surf(x,[0 1], repmat(u(:,nb)-u(:,ne),1,2)'); view(0,90); shading flat;
  u0 = u(:,nb);  tend = t(ne)-t(nb);
%  [ut,tt] = flowdp8('ks01', [N dx sig], u0, 0.0, tend, h);
%  hold on; plot(x,ut,'*-');
end,

if 0,  % Check flow calculation in kszolflow
  f0 = kszolflow(u0, dx, sig);
  [upt,tt] = flowdp8('ks01', [N dx sig], u0, 0.0, del, del);
  disp(f0 - (upt-u0)./del);
end,

if 0,   % Check Jacobian calculation
  % Jacobian from ks01t
  M = eye(N);  uu = [u0 M];  uu = uu(:);  par = [N*(N+1) dx sig N];
  tic; [uut,tt] = flowdp8('ks01t', par, uu, 0.0, tend, h); toc;
  M = reshape(uut(N+1:end),N,N);
  if 0, del = 1e-4;  % Check Jacobian by finite differences
    for ii = 1:N,
      up = u0; up(ii) = up(ii)+del;  um = u0; um(ii) = um(ii)-del;
      [upt,tt] = flowdp8('ks01', [N dx sig], up, 0.0, tend, h);
      [umt,tt] = flowdp8('ks01', [N dx sig], um, 0.0, tend, h);
      figure(3); plot(ii,sum(abs(M(:,ii)-(upt-umt)./(2*del)))./sum(abs(M(:,ii))),'.'); hold on; drawnow;
    end,  end,
end,

if 0,   %%% Locate UPO by Newton-Armijo
%  load ksz66sd4; x = (1:N).*dx;
  NFEVAL = 0;
  for itr = 1:10,
    M = eye(N);  uu = [u0 M];  uu = uu(:);  
    par = [N*(N+1) dx sig N];
    [uut,tt] = flowdp8('ks01t', par, uu, 0.0, tend, h); 
    NFEVAL=NFEVAL+N+1;
    M = reshape(uut(N+1:end),N,N);  ut = uut(1:N);
    f0 = kszolflow(u0, dx, sig);  ft = kszolflow(ut, dx, sig);
    g = u0 - ut;  ng = norm(g);  invDX = inv([M-eye(N) ft; f0' 0]);
    disp(sprintf('%3d  %5d  %8.5f  %12.4e', itr, NFEVAL, tend, ng));
    for ii = 0:18,  % Armijo iteration
      dX = 2.^(-ii).*(invDX*[g; 0]);
      u1 = u0 + dX(1:N);  te1 = tend + dX(N+1);
      if abs(max(u1(:))) < 3.2 & te1 < 2.0*tend & te1 > 0.5*tend,
        [ut,tt] = flowdp8('ks01', [N dx sig], u1, 0.0, te1, h); 
        NFEVAL = NFEVAL+1; g1 = u1 - ut;  
        disp(sprintf('   %2d  %12.4e  %12.4e',ii,norm(dX),norm(g1)));
        if norm(g1) < ng.*(1-2.^(-ii-1)), u0 = u1;  tend = te1; 
          figure(2); plot(x,u0,'r.-',x,ut,'k.-');
%          figure(3); set(hds,'ydata',[get(hds,'ydata') itr+1],...
%            'zdata',[get(hds,'zdata'); u0' - ut']); 
          ii = 123; break; end, end,
    end, if (ii ~= 123) disp('No further improvement'); break; end,
  end,
end,

if 0,  %%% Locate UPO by associated flow  dx/ds = f(x) - x (Euler)
  step = 3e-4;
  for itr = 1:2000,
    g = kszolstb(0.0, [u0; tend]);
    u0 = u0 + step.*g(1:N);  tend = tend + step.*g(N+1);
    disp(sprintf('%4d  %8.3f %11.3e',itr,tend,norm(g(1:N))));
    figure(2); hp1 = plot(x,u0,'r.-',x,u0+g(1:N),'k.-');
%    figure(3); set(hds,'ydata',[get(hds,'ydata') itr+1],...
%      'zdata',[get(hds,'zdata'); -g(1:N)']);  pause;
  end, return;
end,

if 0,  %%% Locate UPO by associated flow  dx/ds = Cg(x) (ode15s)
  global C;  C = eye(N);
%  [v,e] = eig(M);  C(3,3) = -1;  C(5,5) = -1;  C(6,6) = -1;  C = real(v*C*inv(v));
%  [u,s,v] = svd(M);  C = eye(N) + u(:,1:3)*(duposp(3,4)-eye(3))*u(:,1:3)';
%  nsb = 3; isb = 6; cs = zeros(3); cs(1:2,1:2) = duposp(2,isb); cs(3,3) = -1;
%  [u,s,v] = svd(M);  C=eye(N)+v(:,1:nsb)*(cs-eye(nsb))*v(:,1:nsb)';
  global X;  X = x;
  reltol = 1e-7;  abstol = 1e-6;
  options = odeset('abstol',abstol,'reltol',reltol,'outputfcn',@kszolplot);
  [s,y] = ode15s(@kszolstb, [0 20], [u0; tend],options);
end,

if 1,  %%% Locate UPO by associated flow  dx/ds = g(x) from seed
%  load ksz66sd4;
  NFEVAL = 0;
  if 1,
    M = eye(N);  uu = [u0 M];  uu = uu(:);  par = [N*(N+1) dx sig N];
    [uut,tt] = flowdp8('ks01t',par,uu,0.0,tend,h); NFEVAL=NFEVAL+N+1;
    M = reshape(uut(N+1:end),N,N);  ut = uut(1:N);  end,
%    [u,s,v] = svd(M);  [w,t] = schur(M);  end,  return;
  global C;  C = eye(N);  nsb = 2;
  for isb = 0:7,  
%  [v,e] = eig(M);  C(3,3) = -1;  C(5,5) = -1;  C(6,6) = -1;  C = real(v*C*inv(v));
%    cs = zeros(3);  cs(1:2,1:2) = duposp(2,isb); cs(3,3) = -1;
    disp(isb); cs = duposp(nsb,isb);% cs = randn(nsb); [cs,s,v] = svd(cs);
    if 1, [u,s,v] = svd(M);  else,  [u,s] = schur(M);  end,
    C = eye(N) + u(:,1:nsb)*(cs-eye(nsb))*u(:,1:nsb)';
    e = eig(C*(M-eye(N)));  disp(real(e(1:5))');
    global X;  X = x;  reltol = 1e-7;  abstol = 1e-6;
    options = odeset('abstol',abstol,'reltol',reltol,'outputfcn',@kszolplot);
    [s,y] = ode15s(@kszolstb, [0 5], [u0; tend],options); end,
end


if 0,  %%% Locate UPO using Matlab Optimisation Toolbox routines
%  load ksz33sd1v6;  x = (1:N).*dx;
  [u,t] = flowdp8('ks01', [N dx sig], u0, 0.0, (0:.01:2).*tend, h);
  figure(1);  clf; pcolor(x,t,u'); shading flat; hold on; 
  plot(x([1 end]), [tend tend], 'k-');  disp(norm(u(:,1)-u(:,101)));  pause;
  options = optimset('Display','iter','MaxIter',100,'Jacobian','on');
  lb = [-3.2.*ones(size(u0)); 0];  ub = [3.2.*ones(size(u0)); 100];
  [up,gn,gp,exitflag,output] = lsqnonlin(@kszolopt,[u0; tend],lb,ub,options);
end