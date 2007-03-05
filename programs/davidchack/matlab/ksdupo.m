%% Locate close returns in real space
  N = 31;  d = 2*pi/sqrt(0.015);  h = 0.2;  randn('seed',12340014);
  u = 0.1*randn(N,1);
  [tt, uu] = ksedtodd(d, u, h, 1000, 5);
  figure(1); clf;  colormap(jet);
  pcolor(d.*(1:N)'./(2*N+2), tt, uu');  shading flat;
   
  if 1, cr = zeros(size(uu,2));
    for ii = 1:size(uu,2),  cr(ii,:) = sqrt(sum((uu-repmat(uu(:,ii),1,size(uu,2))).^2));  end,
    figure(2); clf; %pcolor(cr); shading flat; hold on; 
    contour(1:201, 1:201, cr,1:3);
  end,


%% Locate close returns in Fourier modes
  N = 63;  nu = (2*pi/54.0)^2;  nu = 0.015;
  h = 0.002;  randn('seed',12340000);
  a0 = randn(N,1).*exp(-0.3*(1:N)'); % a0 = zeros(N,1); a0(1) = -1.0;
  [tt, aa, a0p, ap] = kscvfm(nu, 5.0, a0, h, 5);  na = size(aa,2);
  uu = real(ifft(1i*(2*N+2)*[zeros(1,na); aa; zeros(1,na); -flipud(aa)]));
  uu = uu(2:N+1,:);  xx = pi*(1:N)'/(N+1);
  figure(1); clf; pcolor(xx,tt,uu'); shading flat;
  
  global C;  C = eye(N);

  switch 1,
  case 1,
    cr = zeros(na);
    for ii = 1:size(aa,2),  cr(ii,:) = sqrt(sum((aa-repmat(aa(:,ii),1,na)).^2));
%      cr(ii,:) = cr(ii,[ii:end 1:ii-1]);
      cr(ii,:) = [cr(ii,ii:end) zeros(1,ii-1)];
    end,
    figure(2); clf;
    contourf(tt, tt, cr', 0:1:4);  grid on;  xlabel('t_0');  ylabel('T','rotat',0);
    if 1,  [nb, ne] = ginput(1);  ne = ne + nb;
      nb = round(nb./(tt(2)-tt(1)));  ne = round(ne./(tt(2)-tt(1)));
    else,
%  12340000
%      nb = 215;  ne = 287; % 3: 0 - 0.7294, 1 - 0.660 (almost)
%      nb = 191;  ne = 252; % 3: 0,1 - 0.6826, 2 - 0.6819
%      nb = 123;  ne = 262; % 3: 0 - 1.3936
%      nb = 123;  ne = 203; % 3: 2 - 0.7294
%      nb = 123;  ne = 168; % 3: 0 - 0.5057
%      nb = 275;  ne = 395; % 2: 0 - 1.2646, 3 - 1.2678
%      nb = 158;  ne = 224; % 3: 0,1 - 0.6826
%  12340050
%      nb = 121;  ne = 217; % 3: 0,1 - 0.9205, 7 - 0.9968
%  -1.01
%       nb = 234;  ne = 351; % 3: 0,2 - 1.1845, 1,3 - 1.2677,
%nu = 0.025534
%  12340000
%       nb = 157;  ne = 430;  % 1: 0 - 2.7383,
%  12340001
%       nb =  86;  ne = 176;  % 3: 0 - 0.9961,
%N = 31;  nu = (2*pi/38.5)^2;  h = 0.002;  randn('seed',12340000);
%     nb = 146;  ne = 278; % 3: 0 - 1.4017, 1,2,3 - 1.3821, 8 - 2x0.6826
%N = 63;  nu = (2*pi/54.0)^2;  h = 0.002;  randn('seed',12340000);
%     nb = 54;  ne = 104; % 5: 0 - 0.4939, 
%     nb = 136; ne = 187; % 3: 0,10 - 0.4838, 
%     nb = 256; ne = 272; % 3: 0 - 0.3500,
    end,
    if ne < nb;  ii = nb;  nb = ne;  ne = ii; end,
    a0 = aa(:,nb);  tend = tt(ne)-tt(nb);  disp([nb ne]);% return;
    if 1,  %  Determine stretching directions at the seed
      nud = 3;  nst = 2;
      [ts, as, da] = kscvfmj(nu, tt(nb), aa(:,1), h, 0); 
      [uu, dd, vv] = svd(da);
      disp([find(diag(dd)>1)' nst]);  disp(duposp(nud,nst));
      C = eye(N) + uu(:,1:nud)*(duposp(nud,nst)-eye(nud))*uu(:,1:nud)';
    else
      [tt, aa, da] = kscvfmj(nu, tend, a0, h, 0);
      [uu,dd,vv] = svd(da-eye(N));  C = -vv*uu';
    end,
  case 2,
%    nb = 85;   ne = 124;%  stationary point  
%    nb = 157;  ne = 226;  
%    nb = 180;  ne = 322;
%    nb = 61;   ne = 398;%  tend = 3.363770;
%    nb = 358;   ne = 475;%  tend = 1.184548;
    a0 = aa(:,nb);  tend = tt(ne)-tt(nb);
  case 3,
    load ksfm064upo;  N = 31;  a0 = a0(1:N);
    [tt, aa, da] = kscvfmj(nu, tend, a0, 0.002, 0);
%    [vv,dd] = eig(da);  dd(1,1) = -dd(1,1);  da = real(vv*dd*inv(vv));
    [uu,dd,vv] = svd(da-eye(N));  C = -vv*uu';
  case 4,
%    nb = 83;   ne = 262;
    nb = 308;  ne = 347;
    a0 = aa(:,nb);  tend = tt(ne)-tt(nb);
    [tt, aa, da] = kscvfmj(nu, tend, a0, 0.002, 0);
    [uu,dd,vv] = svd(da-eye(N));  C = -vv*uu';
  case 5,
    load ksfm278seed2;
    [tt, aa, da] = kscvfmj(nu, tend, a0, 0.002, 0);
    [uu,dd,vv] = svd(da-eye(N));  C = -vv*uu';
  case 6,
    load ksfm0920upo;  ip = 1;  a0 = as(:,ip);  tend = ts(end);  tend = 1.5;
    [tt, aa, da] = kscvfmj(nu, tend, a0, h, 0);
    [uu,dd] = eig(da);
    dd(1,1) = -dd(1,1);
%    dd(2,2) = -dd(2,2);
%    a0 = a0 - 1e-3.*uu(:,2);
    das = real(uu*dd*inv(uu));  [uu,dd,vv] = svd(das-eye(N));  C = -vv*uu';
  end
%  return;
  reltol = 1e-7;  abstol = reltol; %[10.*reltol.*exp(-0.3*(1:N)'); 1];
  options = odeset('abstol',abstol,'reltol',reltol,'outputfcn',@ksdupoplot);
  [s,y] = ode15s(@ksfmc5, [0 40], [a0; tend],options);

%% Plot real-space KS-odd solution
  clear; load ksfm064upo;  t = (0:0.005:1)'*tend;
  tic; b = ksfmflowmapj(nu, t, a0, 2e-4); toc;
  d = 2*pi/sqrt(nu);
  y = (0:N)'*pi/N; w = zeros(length(y),length(t));
  for i=1:length(t), w(:,i) = -2*sum(repmat(b(:,i),1,length(y)).*sin(y*(1:N))')'; end,
  u = -4*pi.*w./d;  x = y*d/(2*pi);
  figure(1); clf; pcolor(x,t/nu,u'); shading flat;
  
  u0 = u(2:end-1,1);
  
  tic; [tt, uu] = ksedtodd(d, u0, 0.2, 212, 1); toc;
  figure(2); clf; pcolor(x(2:end-1),tt,uu'); shading flat;


%% test FFT calculation of ksfm coupled term
  d = 7;  a = (1:d)';  k = (1:d);  f = zeros(d,1);  N = 2*(d+1);
  for ik = k,
    f(ik) = sum(a(1:ik-1).*a(ik-1:-1:1)) - 2.*sum(a(1:d-ik).*a(ik+1:d));
  end,
  
  uk = 1i*[0; a; 0; -flipud(a)].*N;
  uu = real(ifft(uk));
  v = -ifft(uu.^2);
  
  aa = [-flipud(a); 0; a];  bb = [0; flipud(aa(2:end))]';
  for ik = k(2:end),
    bb = [bb; [-aa(end-ik+2:end)' 0 flipud(aa(ik+1:end))']];
  end
  
  af = [0; a; 0; -flipud(a)];


%% solve ksfmfj flow with ode15s; identify local minima of norm(f)
  clear all;  L = 50;  d = 62;
  nu = (2*pi/L).^2;  a0 = zeros(d,1);
  randn('seed',12340003);  a0(1:7) = 0.5*randn(7,1);
  [t, a] = ode15s(@ksfmfj, [0 10], a0, [], nu);
  nt = length(t);
  x = (0:d+1)./(d+2)*L;
  u = real(ifft([zeros(1,nt); a(:,1:2:d-1)'+1i*a(:,2:2:d)'; zeros(1,nt); a(:,d-1:-2:1)'-1i*a(:,d:-2:2)'])).*(d+2);
  figure(1);  set(gcf,'pos',[40  80  300  850]);  
  pcolor(x,t,u'); shading flat;  set(gca,'pos',[.13 .04  .82 .92]);
  ng = zeros(nt,1); 
  for ii = 1:nt,
    f = ksfmfj(0, a(ii,:)', nu);  ng(ii) = norm(f);
  end,
  figure(2); clf; set(gcf,'pos',[350 640  1000  300]);
  plot(t, ng,'.-');  set(gca,'pos',[.04 .13 .92 .82]); grid on; hold on;
  im = find(ng(2:end-1)<2 & ng(2:end-1)-ng(1:end-2)<0 & ng(3:end)-ng(2:end-1)>0)+1;
  plot(t(im),ng(im),'ro');  ie = 0; return;
  for jm = 1:1,
    if 1,  % NA approach
      ae = kseqna(nu, a(im(jm),:)');
      ge = ksfmfj(0, ae, nu);
      disp(sprintf('%3d  %10.3e', jm, norm(ge)));
      if norm(ge) < 1e-4,
        ie = ie + 1;
        ue = real(ifft([0; ae(1:2:d-1)+1i*ae(2:2:d); 0; ae(d-1:-2:1)-1i*ae(d:-2:2)]));
        figure(3), plot(x,ue([i0:end 1:i0-1]),'.-'); hold on;
%        [ma,i0] = max(ue);  ve = fft(ue([i0:end 1:i0-1]));
%        figure(4); plot(abs(ve(3)),abs(ve(5)),'.'); hold on;
        figure(2), plot(t(im(jm)),0,'k*');
      end,
    else,  % Polar decomposition approach
      ae = a(im(jm),:)';  ti = [0.35 2.25 100];
      for ii = 1:3,
        [fm,dfm] = ksfmfj(0, ae, nu);
        [u,s,v] = svd(dfm);  C = -v*u';
        [t, a] = ode15s(@ksfmfj, [0 ti(ii)], ae, [], nu, C);
        ae = a(end,:)';
      end,
    end,
  end,

%% Check correspondense between Fortran and Matlab output for ksfmedt
  d = 22;  N = 32;  x = (1:N)'./N*d;  q = 2*pi/d;  u = cos(0.5*q*x);
  v = ifft(u);  a0 = [real(v(2:N/2)) imag(v(2:N/2))]'; a0 = a0(:);
  a0 = zeros(N-2,1); a0(1:6) = 0.2*randn(6,1);
  [tt,aa,a0p,ap] = ksfmedt(d, 500, a0, 0.25, 4);
  figure(1); clf; plot(tt,aa(1:4,:),'.-'); 
  v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
  v = [zeros(1,size(aa,2)); v; zeros(1,size(aa,2)); flipud(conj(v))];
  u = real(fft(v)); figure(2); clf; pcolor(x,tt,u'); shading flat;


%% Locating steady states of ksfm using fsolve
  clear;  load ks22uss3a;  a0 = [a0; zeros(32,1)];
  
  d = 22;  N = 32;  x = (1:N)'./N*d;  h = 0.25;  tpre = 10;  t0 = 50;
%  a0 = zeros(N-2,1);  a0(2) = 1; %randn('seed',12340019);  a0(1:8) = 0.3*randn(8,1);
  load ks22seed001a;
  [tt,aa] = ksfmedt(d, tpre, a0, h, 1);
  nap = zeros(size(aa,2),1);
  for ii = 1:size(aa,2), nap(ii) = norm(ksfm(t0,aa(:,ii),d)); end,
  [map,iap] = min(nap); tpre = tt(iap);
  figure(1); clf; plot(tt,nap,'.-',tpre,map,'ro'); pause;
  [tt,aa] = ksfmedt(d, tpre, a0, h); % aa = a0;
  options = optimset('Display','iter','Jacobian','on','TolFun',1e-9);
  tic; a0 = fsolve(@(a)ksfm(t0,a,d),aa,options); toc;
  [tt,aa] = ksfmedt(d, t0, a0, h, 1);
  nap = zeros(size(aa,2),1);
  for ii = 1:size(aa,2), nap(ii) = norm(ksfm(t0,aa(:,ii),d)); end,
  hold on; plot(tt,nap,'k.-');
  v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
  v = [zeros(1,size(aa,2)); v; zeros(1,size(aa,2)); flipud(conj(v))];
  u = real(fft(v)); figure(2); clf; pcolor(x,tt,u'); shading flat;

%% Locating traveling wave with ksfmtr using fsolve
    clear;  d = 18.5;  N = 32;  h = 0.25;
    x = d.*(1:N)'./N;  np = 1;  tend = 200;
    a0 = zeros(N-2,1); randn('seed',12340003);  a0(1:8) = 0.2*randn(8,1);
    [tt, aa, a0p, ap] = ksfmedt(d, tend, a0, h, np);  na = size(aa,2);
    v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
    v = [zeros(1,size(aa,2)); v; zeros(1,size(aa,2)); flipud(conj(v))];
    u = real(fft(v));  
    figure(1); clf; subplot(1,2,1); pcolor(x,tt,u'); shading flat; 
    k = (2*pi/d).*[0:N/2-1 0 -N/2+1:-1]';
    fa = repmat(k.^2.*(1-k.^2),1,na).*v + ...
         0.5i*repmat(k,1,na).*ifft(u.^2);
    c = (-.5:0.01:.5)';  nc = length(c); ha = zeros(nc,na);
    for ii = 1:nc,
      hc = fa(2:N/2,:) - 1i*c(ii)*repmat(k(2:N/2),1,na).*v(2:N/2,:);
      ha(ii,:) = sqrt(sum(hc.*conj(hc)));
    end,
    subplot(1,2,2);
%    contourf(c,tt,ha',[0:.01:.1]); shading flat; wide(1.1);  grid on;
    [mha,mc] = min(ha);  plot(c(mc),tt,'r.');  grid on; pause;
    ic = 650; 
    ic = 500;
    ic = 48;
    ic = 10;
    options = optimset('Display','iter','Jacobian','off','TolFun',1e-9);
    a0 = fsolve(@(a)ksfmtr(a,d),[aa(:,ic); c(mc(ic))],options); pause;
    for d = 18.5:.5:22;  disp(['L = ' num2str(d)]);
      a0 = fsolve(@(a)ksfmtr(a,d),a0,options);  end,
  
%% Locating TW2
  clear; load ks22seed000a;  c = ph/tend;
  options = optimset('Display','iter','Jacobian','off','MaxFunEvals',20000,'TolFun',1e-9);
  a0 = fsolve(@(a)ksfmtr(a,d),[a0; c],options);
    
%%  Stability of the steady states and traveling wave
  clear; 
%  load ks22uss2b; c = 0;
%   load ks22uss3a; c = 0;
  load ks22utw2a;
  c = -c;  a0(1:2:end) = -a0(1:2:end);
%  a0 = [a0; zeros(32,1)];
  [f, df] = ksfm(0, a0, d);  na = length(a0);
  ik = 2*pi*c/d*(1:na/2)';  h = f;  
  h(1:2:end) = h(1:2:end) + ik.*a0(2:2:end);
  h(2:2:end) = h(2:2:end) - ik.*a0(1:2:end);
  
  dh = df;  
  dh(2:2*na+2:end) = dh(2:2*na+2:end) - ik';
  dh(na+1:2*na+2:end) = dh(na+1:2*na+2:end) + ik';
  [vdh, edh] = eig(dh);  edh = diag(edh);
  [sedh, ie] = sort(real(edh),1,'descend');
  edh = edh(ie);  vdh = vdh(:,ie);  disp(edh(1:10));

  
%% Locate steady states of ksfm using stabilising transformations.
  clear;
  d = 22;  N = 32;  x = (1:N)'./N*d;  h = 0.25;  tpre = 50;  t0 = 50;
  a0 = zeros(N-2,1);  randn('seed',12340021);  a0(1:8) = 0.3*randn(8,1);
  [tt,aa] = ksfmedt(d, tpre, a0, h, 1);
  nap = zeros(size(aa,2),1);
  for ii = 1:size(aa,2), nap(ii) = norm(ksfm(t0,aa(:,ii),d)); end,
  [map,iap] = min(nap); tpre = tt(iap);
  figure(3); clf; plot(tt,nap,'.-',tpre,map,'ro'); pause;
  [tt,aa] = ksfmedt(d, tpre, a0, h);
  [f,df] = ksfm(t0,aa,d); [u,s,v] = svd(df);  C = -v*u';
  options = odeset('Stats','on');
  tic; [tt, a0] = ode15s(@(t,a)ksfm(t,a,d,C),[0 500],aa,options); toc;
  a0 = a0(end,:)';  [tt,aa] = ksfmedt(d, t0, a0, h, 1);
  nap = zeros(size(aa,2),1);
  for ii = 1:size(aa,2), nap(ii) = norm(ksfm(t0,aa(:,ii),d)); end,
  disp(nap(1)); hold on; plot(tt,nap,'k.-');
  v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
  v = [zeros(1,size(aa,2)); v; zeros(1,size(aa,2)); flipud(conj(v))];
  u = real(fft(v)); figure(4); clf; pcolor(x,tt,u'); shading flat;
  

%% Locate UPOs in ksfmedt flow
  clear;
  d = 22;  N = 32;  x = d.*(1:N)'./N;  h = 0.25;  tpre = 300;  t0 = 50;
  a0 = zeros(N-2,1);  randn('seed',12340005);  a0(1:8) = 0.2*randn(8,1);
  [tt, aa, a0p, ap] = ksfmedt(d, tpre, a0, h, 1);  na = size(aa,2);
  v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
  v = [zeros(1,size(aa,2)); v; zeros(1,size(aa,2)); flipud(conj(v))]; u = real(fft(v));
  figure(1);  set(gcf,'pos',[5 35 250 910]); clf; 
  pcolor(x,tt,u'); shading flat;
  switch 1,
  case 1,  %%% Choose close return interactively
    cr = zeros(na);
    for ii = 1:size(aa,2),  cr(ii,:) = sqrt(sum((aa-repmat(aa(:,ii),1,na)).^2));
      cr(ii,:) = [cr(ii,ii:end) zeros(1,ii-1)]; end
    figure(2); set(gcf,'pos',[685 300 710 650]); clf;
    contourf(tt, tt, cr', 0:.4:2);  grid on;  xlabel('t_0');  ylabel('T','rotat',0);
    if 1,  [nb, ne] = ginput(1);  ne = ne + nb;
      nb = round(nb./(tt(2)-tt(1)));  ne = round(ne./(tt(2)-tt(1))); end
    a0 = aa(:,nb);  tend = tt(ne)-tt(nb);  disp([nb ne norm(aa(:,ne)-aa(:,nb))]);% return;
    figure(1); hold on; plot(x([1 end end 1]),tt([nb nb ne ne]),'k-');
  case 2,  %%% Close return after tpre within t0
    tpre = 50;  nb = min(find(tt>tpre));
    na = sqrt(sum((aa(:,nb:end)-repmat(aa(:,nb),1,size(aa(:,nb:end),2))).^2));  dna = diff(na);
    ie = find(dna(1:end-1) < 0 & dna(2:end) > 0);  ne = ie(4)+nb;
    figure(2); clf; plot(tt(nb:end)-tt(nb),na','g.-',tt(ie+nb),na(ie),'ro');  grid on;
    a0 = aa(:,nb);  tend = tt(ne)-tt(nb);  disp([tend na(ne)]);
    figure(1); hold on; plot(x([1 end end 1]),tt([nb nb ne ne]),'k-');
  case 3,  %%% Load seed from a file
    load ks22seed65a;
  end
%  N = 64; x = d.*(1:N)'./N; a0 = [a0; zeros(32,1)];
  switch 2,
  case 1,  %%% associated flow ksfmf1
    C = eye(N-2);  alpha = 10000;  ph = 0;
    options = odeset('reltol',1e-7,'outputfcn',@ksdupoplot);
    [s,y] = ode15s(@(t,a)ksfmf1(t,a,d,h,C,alpha), [0 100], [a0; tend], options);
    sv = input('Save final orbit? (y/n): ','s');
    if sv == 'y', fname = input('File name: ','s'); 
      a0 = y(end,1:end-1)';  tend = y(end,end);
      save(fname,'d','a0','tend'); end
%    a0 = y(end,1:end-2)';  tend = y(end,end-1);  ph = y(end,end);
%    save(fname,'d','a0','tend','ph'); end
  case 2,  %%% associated flow ksfmfm (with T from min |g| )
    global TT A F FT NFEVAL DD
    TT = tend;  C = eye(N-2);  DD = d;
    if 0,  [tt,aa,da] = ksfmjedt(d,tend,a0,h);  e = eig(da);  disp(e(1:5))
      switch 2,
      case 1, [u,C,v] = svd(da);  C = eye(N-2) - 2*u(:,1)*u(:,1)';
      case 2, [u,C,v] = svd(da-eye(N-2));  C = -v*u';
      end
    end
    options = odeset('abstol',1e-9,'reltol',1e-10,'outputfcn',@ksdupoplot3);
    [s,y] = ode15s(@(t,a)ksfmfl(t,a,d,h,C), [0 100], a0, options);
    sv = input('Save final orbit? (y/n): ','s');
    if sv == 'y', fname = input('File name: ','s');
      a0 = y(end,:)'; tend = TT(end); save(fname,'d','h','a0','tend'); end
  end

%% Locate UPOs of ksfmedt with constraint a(1) = 0 - doesn't work
  clear;
  d = 22; N = 32;  x = d.*(1:N)'./N;  h = 0.25;  tpre = 200;  t0 = 100;
  a0 = zeros(N-2,1);  randn('seed',12340000);  a0(2:9) = 0.1*randn(8,1);
  [tt, aa] = ksfmedt(d, tpre, a0, h, 1);
  figure(1); clf; plot(tt,aa(1:2,:),'.-');  grid on;
  i0 = find(aa(1,1:end-1) < 0 & aa(1,2:end) > 0);  ii0 = 3;
  a0 = aa(:,i0(ii0));  a0(1) = 0; 
  [tt, aa] = ksfmedt(d, t0, a0, h, 1);
  na = sqrt(sum((aa-repmat(aa(:,1),1,size(aa,2))).^2));  dna = diff(na);
  ie = find(dna(1:end-1) < 0 & dna(2:end) > 0);  iie = 2;
  figure(2); clf; plot(tt,na','g.-',tt(ie+1),na(ie+1),'ro');  grid on;
  tend = tt(ie(iie)+1); disp([tend na(ie(iie)+1)]);
  
  C = eye(N-3);  alpha = 5;
  options = odeset('reltol',1e-7,'outputfcn',@ksdupoplot);
  [s,y] = ode15s(@(t,a)ksfmf2(t,a,d,h,C,alpha), [0 1000], [a0(2:end); tend], options);


%% Locate UPOs of ksfmedt with fsolve (from close return seed)
  clear; load ks22seed66c;  h = 0.25;
  global NFEVAL hp1 hp2;  NFEVAL = 0;
  [tt, aa] = ksfmedt(d, tend, a0, h, 2);
  figure(3); clf; hp1 = plot(aa(1,:),aa(2,:),'.-'); grid on; hold on;
  plot(aa(1,1),aa(2,1),'ro',aa(1,end),aa(2,end),'ko');
  figure(5); clf;  N = size(aa,1)+2;  x = d.*(1:N)'./N;  
  v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
  v = [zeros(1,size(aa,2)); v; zeros(1,size(aa,2)); flipud(conj(v))];
  u = real(fft(v)); hp2 = pcolor(x,tt,u'); shading flat; axis tight;
  options = optimset('Display','iter','MaxFunEvals',15000); %,'Jacobian','on','TolFun',1e-9);
  upo = fsolve(@(a)ksfmf3(a,d,h),[a0;tend],options);

%% Use UQOs as seeds
  clear;  d = 22;  N = 32;  x = d.*(1:N)'./N;  h = 0.25;
  load ks22uqo042a;  te = 200;
  [tt, aa, da] = ksfmjedt(d, tend, a0, h, 2);
  k = (-2i*pi/d).*(1:N/2-1)';  ek = exp(k.*ph);
  cda = repmat(ek,1,N-2).*(da(1:2:end,:)+1i*da(2:2:end,:));
  da(1:2:end,:) = real(cda);  da(2:2:end,:) = imag(cda); 
  [u,C,v] = svd(da-eye(N-2));  C = -v*u';
  eda = eig(C*(da-eye(N-2)));  
  [tt2, aa2, da2] = ksfmjedt(d, 2*tend, a0, h, 4);
  ek2 = exp(2.*k.*ph);
  cda = repmat(ek2,1,N-2).*(da2(1:2:end,:)+1i*da2(2:2:end,:));
  da2(1:2:end,:) = real(cda);  da2(2:2:end,:) = imag(cda);
  [v,e] = eig(da);  s = eye(N-2);  s(1,1) = -1;
  [u,C,v] = svd(v*s*e*inv(v)-eye(N-2));  C = -v*u';
  eda2 = eig(C*(da2-eye(N-2)));

%%
  v = ek.*(aa(1:2:end,end)+1i*aa(2:2:end,end));
  ap = a0;  ap(1:2:end) = real(v);  ap(2:2:end) = imag(v);
  figure(3); clf; 
  plot([aa(1,:) ap(1)],[aa(2,:) ap(2)],'.-',aa(1,1),aa(2,1),'ko');  grid on;
  [tt, aa] = ksfmedt(d, te, a0, h, 2);  
  figure(2);  clf;  phase = (-0.5:0.01:0.49)'*d;
  v = aa(1:2:end,:)+1i*aa(2:2:end,:);
  fac = exp(k*phase');  cr = zeros(size(fac,2),size(v,2));
  for ii = 1:size(fac,2),
    cr(ii,:) = sqrt(sum(abs(v.*repmat(fac(:,ii),1,size(v,2)) - ...
                    repmat(v(:,1),1,size(v,2))).^2)); end
  contourf(tt,phase,cr,[0:.2:1]); shading flat;
  [tend, ph] = ginput(1);  ek = exp(k.*ph);
  cda = repmat(ek,1,N-2).*(da(1:2:end,:)+1i*da(2:2:end,:));
  da(1:2:end,:) = real(cda); da(2:2:end,:) = imag(cda); 
  eda = eig(da); [me,ie] = sort(abs(eda),'descend'); disp(eda(ie(1:5))');
  [u,C,v] = svd(da-eye(N-2));  C = -v*u';

%% Convert real-space initial condition to fourier modes
% function a0 = rs2fm(u0)
  clear; d = 13;  N = 32;  sig = 0.2;  np = 4;  h = 0.25;
  x = d.*(1:N)'./N;  u0 = (x-d/2).*exp(-(x-d/2).^2/(2*(2*pi*sig).^2));  u0(end) = 0;
  v0 = ifft(u0);  a0 = zeros(N-2,1);  
  a0(1:2:end) = real(v0(2:N/2));  a0(2:2:end) = imag(v0(2:N/2));
  a0 = zeros(N-2,1); a0(1:8) = 0.1*randn(8,1);
  [tt, aa] = ksfmedt(d, 100, a0, h, np);  na = size(aa,2);
  v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
  v = [zeros(1,size(aa,2)); v; zeros(1,size(aa,2)); flipud(conj(v))];
  u = real(fft(v)); 
  figure(1); clf; pcolor(x,tt,u'); shading flat; tall(1.1); hold on;
  if 0,  figure(2); clf; 
    for ii = 1:size(u,2),  plot(x,u(:,ii),'.-'); 
      title(['Time: ' num2str(h*np*(ii-1))]); grid on; pause; end, end

%% Locate UQOs of ksfmedt
  clear;  load kse22orbits;  d = 22;  N = 32;  h = 0.25;  
  switch 1,  %%% Switch for seed selection
  case 1,  %%% Find close returns with a shift
    x = d.*(1:N)'./N;  np = 6;  tpre = 500;  tb = 200;  te = 200;  iav = 0;
    if 0,
      a0 = zeros(N-2,1); randn('seed',22000105);  a0(1:10) = 0.2*randn(10,1);
    else
%      a0 = 0.1*ones(N-2,1); end,
      load ks22utw1a; a0(1:2:end) = -a0(1:2:end);
%      load ks22uss2b; a0(3) = a0(3)+1e-2;
%      load ks22seed062b;
%      load ks22uqo170b;
%      load ks50rpo070a;
%      load ks22uss2a;
%      load ks22seed2w3w; a0 = a1; clear a1;
    end
    [tt, aa, a0p, ap] = ksfmedt(d, tpre, a0, h, np);  na = size(aa,2);
    v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
    v = [zeros(1,size(aa,2)); v; zeros(1,size(aa,2)); flipud(conj(v))];
    u = real(fft(v));
    figure(1); set(gcf,'pos',[5 35 250 910]); clf; 
    pcolor(x,tt,u'); shading flat; tall(1.1); hold on;
    hp1 = plot(x([1 end]),tt([1 1]),'k-');  
    figure(2); set(gcf,'pos',[263 344 1135 600]); clf;
    for ib = 0:1:200,
      nb = min(round(tb/(h*np))+ib, length(tt));
      ne = min(round(te/(h*np))+nb+1, length(tt));
      set(hp1,'xdata',[x([1 end]);NaN;x([end 1])],'ydata',[tt([nb nb]) NaN tt([ne ne])]);
      a0 = aa(:,nb);  phase = (-0.51:0.01:0.51)'*d;
      v = aa(1:2:end,nb-iav:ne+iav)+1i*aa(2:2:end,nb-iav:ne+iav);  k = (1:N/2-1)';
      fac = exp((-2i*pi./d).*k*phase');  cr = zeros(size(fac,2),size(v,2)-2*iav);
      for ii = 1:size(fac,2), for iia = -iav:2:iav,
          cr(ii,:) = cr(ii,:) + ...
            sqrt(sum(abs(v(:,iia+iav+1:end+iia-iav).*repmat(fac(:,ii),1,size(v,2)-2*iav) - ...
            repmat(v(:,iia+iav+1),1,size(v,2)-2*iav)).^2))...
          ./sqrt(sum(abs(repmat(v(:,iia+iav+1),1,size(v,2)-2*iav)).^2))...
          ./(iav+1); end, end
%          ./sqrt(sum(abs(repmat(v(:,iia+iav+1),1,size(v,2)-2*iav)).^2))...
      contourf(tt(nb:ne)-tt(nb),phase,cr,[0:.1:1]); shading flat;
%      contourf(phase,tt(nb:ne)-tt(nb),cr',[0:.1:1]); shading flat;
      hold on; plot([1;1]*[rpo.T],[1;-1]*[rpo.d],'w.','markersize',16); 
      set(gca,'pos',[0.06 0.08 0.90 0.86]); hold off;
      [tend, ph, but] = ginput(1); disp(sprintf('%4d %2d %6.2f',ib,but,tt(nb))); 
      if but == 32, hold on; plot(tend,ph,'k.','markersize',8); hold off; break; end, end
  case 2,  %%% Load seed from a file
%    load ks50rpo026b; %tend = 167.8; %ph = 5.4988;
    load ks22uqo183b;
%    load ks22upo138a;
%    load ks22uss0a;
%    load ks22seed151a; %tend = 125.2;  ph = 5.26;
%    load ks22seed177a;
  end;
  for sw = 2, %[2,1,3],
  switch sw, %%% Switch for method selection
  case 1,   %%% Find UQOs by associated flow
    alpha = [1e0 1e0];
    for ic = 1:20,
      switch 2,
      case 0,  C = eye(N-2);
      case 1, format short g;
        [g, jac] = ksfmf5(0,[a0;tend;ph],d,h,eye(N),alpha);
        disp(['|g|=' num2str(norm(g(1:end-2)))]);  de = eig(jac(1:end-2,1:end-2)+eye(N-2)); 
        [se,ie]=sort(abs(de),'descend');  disp(de(ie(1:10)));
        [u,C,v] = svd(jac(1:end-2,1:end-2)); C = blkdiag(-v*u',eye(2)); format short;
      case 2, format short g;
        [g, jac] = ksfmf5(0,[a0;tend;ph],d,h,eye(N),alpha);
        disp(['|g|=' num2str(norm(g(1:end-2)))]);  de = eig(jac(1:end-2,1:end-2)+eye(N-2)); 
        [se,ie]=sort(abs(de),'descend'); disp(de(ie(1:10))); [u,C,v] = svd(jac);
        if 1, atol = 1e-5./max(diag(C)); rtol = 10*atol;
          disp(['atol = ' num2str(atol) '  rtol = ' num2str(rtol)]);
        else atol = 1e-13;  rtol = 1e-12; end
        C = -v*u'; format short; end,
      options = odeset('abstol',atol,'reltol',rtol,'outputfcn',@ksdupoplot4);
      [s,y] = ode15s(@(t,a)ksfmf5(t,a,d,h,C,alpha), [0 100000], [a0; tend; ph], options);
      a0 = y(end-1,1:end-2)'; tend = y(end-1,end-1), ph = y(end-1,end),
      but = input('Stop? ','s');  if but == 'y', break; end, end
    figure(2); hold on;  plot(tend,ph,'ko','markersize',7); hold off;
  case 2,   %%% Find UQOs using fsolve
    C = eye(N-2);  alpha = [1e0 1e0];  te = 0;
    options = optimset('Display','iter','outputfcn',@ksoutfun2,'MaxIter',5000, ...
      'MaxFunEvals',100000,'Jacobian','on','NonlEqnAlgorithm','lm','TolX',1e-12,'TolFun',1e-12);
    [upo,fval,eflag] = fsolve(@(a)ksfmf5new(te,a,d,h,C,alpha),[a0;tend;ph],options);
    a0 = upo(1:end-2); tend = upo(end-1), ph = upo(end),
    figure(2); hold on; plot(tend,ph,'ko','markersize',7); hold off;
  case 3,   %%% Use multiple shooting
    ni = [round(tend/12) 0];  ni(2) = ceil(tend./(h.*(ni(1)+1)));
    disp(['Using multiple shooting: ' num2str(ni)]);
    [tt, aa] = ksfmedt(L, tend, a0, h, ni(2)); ni(1) = ni(1)-1;
    ai = aa(:,1:end-2);  ai = [ai(:); tend-ni(1)*ni(2)*h; ph];
    gi = ksfmfns(ai, L, h, ni);  disp(['|gms| = ' num2str(norm(gi))]);
    options = optimset('Display','iter', 'outputfcn',@ksoutfun2,'MaxIter',2000, ...
      'MaxFunEvals',500000,'Jacobian','off','NonlEqnAlgorithm','lm','TolX',1e-12,'TolFun',1e-12);
    upo = fsolve(@(ai)ksfmfns(ai,L,h,ni),ai,options);
    a0 = upo(1:(length(upo)-2)/(ni(1)+1));  
    tend = upo(end-1)+ni(1)*ni(2)*h;  ph = upo(end);
    figure(2); hold on; plot(tend,ph,'ko','markersize',7); hold off;
  case 4,   %%% Find UQOs using lsqnonlin
    C = eye(N-2);  alpha = [1 1];
    options = optimset('Display','on','outputfcn',@ksoutfun,'MaxIter',2000,'Jacobian','on');
    lb = [repmat(-5,N-2,1);0;-d];  ub = [repmat(5,N-2,1);300;d];
    upo = lsqnonlin(@(a)ksfmf5new(te,a,d,h,C,alpha),[a0;tend;ph],lb,ub,options);
  end, end


%% List all located UPOs and UQOs
  clear;  
  lst = dir('ks22uqo183*.mat');  nst = length(lst);
%  lst = dir('ks50rpo00*.mat');  nst = length(lst);
  for ist = 1:nst,  load(lst(ist).name);  N = length(a0)+2;
    g = ksfmf5new(tend, [a0; tend; ph], d, h, eye(N-2), [1 1]);
    disp(sprintf('    ''%15s'';%% %10.6f %10.6f %12.3e',lst(ist).name, tend, ph, norm(g(1:end-2))));  end

%% Plot all located UPOs and UQOs
  switch 22,
  case 22,  
    clear;  name = [
%    'ks22uqo008a.mat';%   7.889874   5.812044   9.647e-016  <-- RE1

    'ks22uqo016a.mat';%  16.316170   2.863787   1.726e-014
    'ks22uqo020a.mat';%  20.506740   0.000000   1.844e-007
    'ks22uqo033b.mat';%  32.806180  10.956763   9.982e-014
    'ks22uqo034a.mat';%  33.501034   4.045869   2.419e-013
    'ks22uqo035a.mat';%  34.639120   9.599483   1.825e-014
    'ks22uqo036a.mat';%  36.016825   3.911040   3.362e-011

    'ks22uqo036d.mat';%  36.205721   3.813489   1.018e-013
    'ks22uqo040b.mat';%  39.770736  -1.610607   1.812e-008
    'ks22uqo040d.mat';%  40.145637  -3.699441   8.866e-013
    'ks22uqo042a.mat';%  41.507477   6.479327   2.200e-008
    'ks22uqo045a.mat';%  44.709656 -10.605822   1.086e-013
    'ks22uqo046b.mat';%  46.052353   7.624041   2.406e-013

    'ks22uqo047c.mat';%  46.502291  -7.758333   1.318e-013
    'ks22uqo047a.mat';%  47.323866   0.339464   8.331e-012
    'ks22uqo048a.mat';%  47.639322   5.675918   5.407e-010
    'ks22uqo054a.mat';%  53.866635  -2.328494   2.062e-013
    'ks22uqo056a.mat';%  55.599640  -5.247831   4.996e-014
    'ks22uqo060b.mat';%  59.669128   3.979569   4.129e-013

    'ks22uqo060a.mat';%  59.892457  -5.442530   2.973e-012
    'ks22uqo065b.mat';%  64.999343   4.719593   5.821e-014
    'ks22uqo065a.mat';%  65.193986  -6.373938   1.777e-012
    'ks22uqo066b.mat';%  65.540700  -0.582016   3.726e-009
%    'ks22uqo066a.mat';%  65.612385   0.086473   1.671e-014  (33b x 2)
    'ks22uqo067c.mat';%  66.779339   0.000003   3.074e-014
%    'ks22uqo067d.mat';%  67.002067   8.091738   4.233e-014  (34a x 2)
    'ks22uqo068a.mat';%  67.939003   5.566339   6.073e-014

    'ks22uqo068e.mat';%  68.223923   1.706300   2.994e-008
    
    'ks22uqo069a.mat';%  68.537305   4.925039   4.496e-014
    'ks22uqo069b.mat';%  69.065162  -9.700812   3.349e-013
    'ks22uqo070a.mat';%  70.190722   0.721829   8.305e-012
    'ks22uqo072a.mat';%  71.678083   5.502827   1.717e-014
%    'ks22uqo072b.mat';%  72.411828  -7.626821   2.400e-012  (36d x 2)
    'ks22uqo074a.mat';%  73.980562   0.678530   8.188e-013  
    'ks22uqo074c.mat';%  74.380749   4.434083   6.562e-009
    
    'ks22uqo075b.mat';%  75.327068  -3.583450   3.522e-013
    'ks22uqo075a.mat';%  75.483577  -6.235556   2.078e-009  
    'ks22uqo077a.mat';%  76.556365  -3.383947   3.993e-009
    'ks22uqo077c.mat';%  76.604718  -1.677121   1.458e-014
    'ks22uqo078c.mat';%  77.831846 -10.791235   4.383e-014
    'ks22uqo078a.mat';%  78.194219   4.842480   3.825e-008
    
    'ks22uqo079b.mat';%  78.814091   6.710732   5.464e-014
    'ks22uqo079a.mat';%  79.508303   8.872265   4.198e-009
    'ks22uqo080b.mat';%  80.000539   1.802967   4.775e-013
    'ks22uqo080a.mat';%  80.253902  -8.305024   1.442e-013
    'ks22uqo081b.mat';%  80.568032   3.753738   1.477e-013
    'ks22uqo081a.mat';%  80.868624   1.049921   8.927e-009
    
    'ks22uqo084a.mat';%  84.086527 -10.061069   2.248e-013
    'ks22uqo084b.mat';%  84.439682  -5.512903   3.416e-013
    'ks22uqo085a.mat';%  84.989461   4.275548   3.647e-010
    'ks22uqo086a.mat';%  85.744395   4.680305   4.288e-013
    'ks22uqo086c.mat';%  85.840717  -5.502822   9.167e-015
    'ks22uqo086d.mat';%  86.068044  -5.126499   2.367e-014
    
    'ks22uqo086b.mat';%  86.230820 -10.642809   6.009e-014
    'ks22uqo087b.mat';%  86.880781   5.530195   1.565e-012
    'ks22uqo087a.mat';%  87.236323   0.000001   5.367e-014
    'ks22uqo088a.mat';%  88.475386 -10.078385   1.277e-013
    'ks22uqo089c.mat';%  88.616419  -6.209569   9.673e-014
    'ks22uqo089b.mat';%  89.397624  -0.472671   1.560e-012
    
    'ks22uqo089a.mat';%  89.419338  -0.788369   1.031e-007
    'ks22uqo091a.mat';%  91.354619   0.083919   1.993e-008
    'ks22uqo094a.mat';%  94.045042  -6.395471   1.020e-010
    'ks22uqo095a.mat';%  95.252982  -0.000001   1.432e-013
    'ks22uqo096b.mat';%  96.053538  -3.774012   5.968e-013
    'ks22uqo096a.mat';%  96.421035 -10.911410   9.371e-009
    
    'ks22uqo098a.mat';%  98.254857  10.847250   4.403e-012
    'ks22uqo098b.mat';%  98.446990  -1.814935   6.338e-012
    'ks22uqo100a.mat';%  99.612756 -10.969248   8.981e-014
    'ks22uqo100b.mat';%  99.877490  -0.519880   3.599e-012
    'ks22uqo100d.mat';%  99.985049   5.972741   2.096e-012
    'ks22uqo100c.mat';% 100.447577  -5.435634   1.412e-013
    
    'ks22uqo101a.mat';% 101.181669  -2.794991   8.380e-014
    'ks22uqo102a.mat';% 102.789193  -9.206446   1.026e-008
    'ks22uqo103a.mat';% 103.013668   1.006582   3.270e-013
    'ks22uqo103b.mat';% 103.066535   6.743595   2.427e-013
    'ks22uqo103c.mat';% 103.160875   3.911896   1.539e-012
    'ks22uqo105c.mat';% 104.537756   4.837363   1.495e-012
    
    'ks22uqo105a.mat';% 105.426220   0.212326   3.732e-012
    'ks22uqo105b.mat';% 105.489600   7.224473   3.984e-013
    'ks22uqo107a.mat';% 107.269635  10.336297   3.588e-010
    'ks22uqo108b.mat';% 107.681157  -4.491332   1.956e-012
    'ks22uqo108a.mat';% 108.078928 -10.270938   1.988e-011
    'ks22uqo109b.mat';% 109.115133  -6.089034   4.053e-008
    
    'ks22uqo109d.mat';% 109.288714   0.184272   3.416e-010
    'ks22uqo110a.mat';% 109.961306  -5.711897   1.547e-007 
    'ks22uqo110b.mat';% 110.253208   0.790874   1.053e-009
    'ks22uqo111b.mat';% 110.881760  -3.069197   2.070e-013
    'ks22uqo111a.mat';% 111.256292   7.788527   7.935e-013
    'ks22uqo112a.mat';% 111.535813  -0.219776   1.327e-013
    
    'ks22uqo112b.mat';% 112.031992   5.448170   1.488e-013
    'ks22uqo113b.mat';% 112.684863  -2.182842   1.782e-013
    'ks22uqo113c.mat';% 112.870371   7.661806   5.416e-013
    'ks22uqo113a.mat';% 113.664192  -9.847374   8.846e-014
    'ks22uqo115a.mat';% 114.769360   4.912638   4.734e-013
    'ks22uqo115b.mat';% 115.149392  -0.453798   1.362e-013

    'ks22uqo117a.mat';% 116.905525   6.640941   4.431e-010
    'ks22uqo117b.mat';% 117.437473   9.914628   7.289e-013
    'ks22uqo118a.mat';% 117.895472  -7.474551   5.115e-009
    'ks22uqo120a.mat';% 119.556928   8.312556   4.791e-010
    'ks22uqo120b.mat';% 119.947182   6.921632   9.420e-012
    'ks22uqo120c.mat';% 120.064292   7.352324   2.663e-013
    
    'ks22uqo120d.mat';% 120.279488   6.445854   1.266e-013
    'ks22uqo122a.mat';% 121.716005  10.740138   4.017e-010
    'ks22uqo123a.mat';% 122.823942  10.416812   4.744e-013
    'ks22uqo123b.mat';% 123.136374   1.651543   1.480e-012
    'ks22uqo125a.mat';% 125.021554   9.062346   1.270e-009
    'ks22uqo126a.mat';% 126.193329  -3.992781   2.909e-009
    
    'ks22uqo127c.mat';% 126.756021  11.536072   1.929e-012
    'ks22uqo127a.mat';% 127.126195   7.157910   9.227e-013
    'ks22uqo127b.mat';% 127.290946   2.031982   7.409e-012
    'ks22uqo128b.mat';% 127.526866   0.826251   1.674e-012
    'ks22uqo128a.mat';% 127.626882   0.662288   2.362e-012
    'ks22uqo128c.mat';% 128.371655   6.586471   1.420e-011
    
    'ks22uqo129a.mat';% 129.485734  -2.139248   1.103e-011
    'ks22uqo130a.mat';% 130.438430  -5.047473   3.046e-013
    'ks22uqo132b.mat';% 131.981245   1.342006   1.397e-012    

    'ks22upo133a.mat';% 132.638412  -0.000001   3.223e-014

    'ks22uqo133a.mat';% 132.696432  -0.015114   7.568e-013
    'ks22uqo133b.mat';% 133.030048   5.310700   3.773e-012
    'ks22uqo136a.mat';% 135.614655 -10.092283   1.998e-012
    'ks22uqo136c.mat';% 135.715849   2.047110   4.442e-013
    'ks22uqo136b.mat';% 136.065427   3.834442   2.098e-013

    'ks22upo138a.mat';% 138.163981  -0.000005   2.113e-013

    'ks22uqo139b.mat';% 138.903098  -1.530682   4.172e-012
    'ks22uqo139a.mat';% 139.219781   5.790933   4.777e-012
    'ks22uqo139c.mat';% 139.256549   5.389937   9.505e-014
    'ks22uqo140a.mat';% 140.031772   5.136071   5.650e-009
    'ks22uqo142b.mat';% 141.609373  -4.853878   9.659e-013
    
    'ks22uqo142a.mat';% 142.178064  -3.341685   1.220e-009
    'ks22uqo144b.mat';% 143.869603  -4.447408   1.558e-013
    'ks22uqo144c.mat';% 144.038979  -5.515201   3.655e-012
    'ks22uqo144a.mat';% 144.306458  -4.708727   2.288e-013
    'ks22uqo146a.mat';% 145.647372 -10.635471   6.209e-012
    'ks22uqo146b.mat';% 146.340000  -5.893415   3.383e-013

    'ks22uqo148a.mat';% 147.874385  -4.327200   1.776e-009
    'ks22uqo149a.mat';% 149.084398   8.408786   3.700e-012
    'ks22uqo150a.mat';% 150.141052   5.674093   8.729e-012
    'ks22uqo151a.mat';% 151.446067   5.480539   2.275e-012
    'ks22uqo152a.mat';% 152.257156   0.419885   2.179e-012
    'ks22uqo156a.mat';% 156.257970   2.925257   5.687e-012

    'ks22uqo157a.mat';% 157.042870   2.961886   8.743e-012
    'ks22uqo158a.mat';% 158.062421   0.596793   2.928e-010
    'ks22uqo163b.mat';% 163.229756  -3.229459   1.343e-012
    'ks22uqo163a.mat';% 163.423996  -5.979334   2.080e-011
    'ks22uqo164a.mat';% 164.117594 -10.872756   5.843e-012
    'ks22uqo167a.mat';% 166.689318 -10.090257   2.105e-012

    'ks22uqo168a.mat';% 167.549441  -0.746615   4.408e-012    
    'ks22uqo170b.mat';% 169.518906   5.916267   1.308e-010
    'ks22uqo170c.mat';% 170.008790   0.510044   7.809e-012
    'ks22uqo170a.mat';% 170.177047  -8.338197   5.753e-011
    'ks22uqo171a.mat';% 170.936876   6.562272   1.063e-009
    'ks22uqo174a.mat';% 173.888171  -6.655034   2.042e-009

    'ks22uqo175a.mat';% 175.381385   5.432786   5.464e-013
    'ks22uqo181a.mat';% 181.401640  -5.800474   2.933e-012
    'ks22uqo183b.mat';% 182.511632   0.021412   6.808e-013
    'ks22uqo183a.mat';% 182.590494   3.145397   2.067e-011
    'ks22uqo184a.mat';% 184.444306   9.457250   5.639e-012
    'ks22uqo188a.mat';% 188.498020   4.526301   6.670e-013

    'ks22upo194a.mat';% 194.274475   0.000005   2.198e-006

    'ks22uqo196a.mat';% 195.830366   6.488804   4.627e-011
   ];  save ks22names name;  nax = 6;
  case 50,
    clear;  name = [
    'ks50rpo009a.mat';%   8.998497  19.623471   1.700e-015
    'ks50rpo011a.mat';%  11.092722  24.024793   3.651e-015
    'ks50rpo013a.mat';%  13.218836  -4.807292   4.609e-014
    'ks50rpo016b.mat';%  15.557797  20.957970   6.501e-015
    'ks50rpo016a.mat';%  15.669563   9.486030   2.502e-014
    'ks50rpo018a.mat';%  18.182769 -19.657325   1.512e-014
    'ks50rpo020a.mat';%  20.065607 -20.198905   1.677e-014
    'ks50rpo021a.mat';%  20.819836  -1.735221   2.666e-014
    'ks50rpo024a.mat';%  24.104919 -10.811604   6.901e-015
    'ks50rpo026b.mat';%  26.236886 -13.718069   1.676e-014
    'ks50rpo026a.mat';%  26.464834 -19.193295   1.335e-014
    'ks50rpo027a.mat';%  27.388358  10.786236   1.206e-009
    'ks50rpo034b.mat';%  33.594548  -3.156112   1.016e-012
    'ks50rpo034a.mat';%  34.449770   6.355855   4.935e-014
    'ks50rpo039a.mat';%  39.076672   6.512421   1.800e-013
    'ks50rpo041a.mat';%  41.067064   0.323054   5.223e-014
    'ks50rpo047a.mat';%  47.059378 -18.613127   5.588e-014
    'ks50rpo070a.mat';%  70.288203 -12.358193   6.576e-012
    ];  save ks50names name;  nax = 3;
  end
  nst = size(name,1);  np = 2;  te = 200;%  nst = 24;
%  for ist = 1:nst,  ost = 1:nst;
%  nax = 5;  for ist = 1:nax, ost = 2:6, % ost = [3 22 50 57 95];
  for ist = 1:nax,  ost = [73:78];
    switch 1,
    case 1,
      ifig = floor(ist/nax)+1;  iorb = mod(ist-1,nax)+1;
      if iorb == 1,  figure(ifig); clf; orient landscape; 
%        set(gcf,'pos',[35 200 1100 740]);
        set(gcf,'pos',[35 200 1000 640]);
        hax = subplots(1,nax,[0.03 0.03 0.07 0.1],[0.015 0 0 0]); end
      load(name(ost(ist),:));  N = length(a0)+2;  x = d.*(-N:N)'./(2*N);
      ne = ceil(te./tend);
      [tti, aai] = ksfmedt(d, tend, a0, h, np);
      ek = exp((2i*pi/d).*ph.*(1:N/2-1)');
      tt = tti(1:end-1);  aa = aai(:,1:end-1);
      for ie = 1:ne-1,
        vi = (aai(1:2:end,:)+1i*aai(2:2:end,:)).*repmat(ek,1,size(aai,2));
        aai(1:2:end,:) = real(vi);  aai(2:2:end,:) = imag(vi);
        aa = [aa aai(:,1:end-1)];  tt = [tt tti(1:end-1)+ie*tend]; end
      v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
      v = [zeros(1,size(aa,2)); v; zeros(N+1,size(aa,2)); flipud(conj(v))];
      u = real(fft(v));  u = [u; u(1,:)];
      axes(hax(iorb)); pcolor(x,tt,u'); shading flat; caxis([-3 3]);  hold on;
      plot(x([1 end])*ones(1,ne),[tend;tend]*(1:ne),'w-');
      plot(mod([ph;ph]*(1:ne-1),d)-d/2,[(1:ne-1);(2:ne)]*tend,'w-');
      title([sprintf('%3d: T=%6.2f ',ost(ist),tend) '\Delta=' sprintf('%7.3f',ph)]);
  %    title(sprintf('d/T = %8.4f',ph./tend));
      set(gca,'ticklength',0.7*get(gca,'ticklength'),'tickdir','out',...
        'ylim',[0 te],'ytick',20:20:te,'box','off');  xlabel('x');
      if iorb > 1, set(gca,'yticklabel',[]); end
      if iorb == nax,
        hax(nax+1) = axes('pos',[.06 .9 .9 .07]);
        tit1 = 'Kuramoto-Sivashinsky: u_t = -uu_x - u_{xx} - u_{xxxx};  ';
        tit2 = ['x \in [-L/2, L/2];    BC: u(x+L,t) = u(x,t);' sprintf('   L = %4.1f;    ',d)];
        tit3 = 'Solutions of the form: u(x+\Delta,T) = u(x,0)';
        text(0.1,.9,[tit1 tit2 tit3]);  set(hax(nax+1),'visi','off');
        axes(hax(1)); ylabel('t','rotat',0); end
    case 2,  % Plot FMs of RPOs
      load(name(ost(ist),:));  N = length(a0)+2;  x = d.*(-N:N)'./(2*N);
      [tt, aa] = ksfmedt(d, tend, a0, h, 1);  
      v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
      ek = exp((-2i*pi/d).*ph./tend*(1:N/2-1)'*tt);  vr = ek.*v;
      vv = [zeros(1,size(aa,2)); v; zeros(N+1,size(aa,2)); flipud(conj(v))];
      u = real(fft(vv));  u = [u; u(1,:)];
      fig1 = figure('PaperPosition',[0.6345 6.345 20.3 15.23],...
        'PaperSize',[20.98 29.68],'Position',[5 200 900 740]);
      ax(1) = axes('Position',[0.04 0.08 0.13 0.86],'Parent',fig1);
      ax(2) = axes('Position',[0.24 0.54 0.33 0.39],'Parent',fig1);
      ax(3) = axes('Position',[0.64 0.54 0.33 0.39],'Parent',fig1);
      ax(4) = axes('Position',[0.24 0.08 0.33 0.39],'Parent',fig1);
      ax(5) = axes('Position',[0.64 0.08 0.33 0.39],'Parent',fig1);
      axes(ax(1));
      pcolor(x,tt,u'); shading flat; hold on;
      title(sprintf('T = %6.2f  d = %9.6f',tend,ph));  
      xlabel('x'); ylabel('Time');
      for iax = 1:4,
        axes(ax(iax+1));
        plot(real(vr(iax,:)),imag(vr(iax,:)),'.-');  hold on;
        plot(real(vr(iax,1)),imag(vr(iax,1)),'k.','markersize',15);  hold on;      
        plot(real(vr(iax,3)),imag(vr(iax,3)),'r.','markersize',15);  hold on;      
        plot(real(vr(iax,20:20:end)),imag(vr(iax,20:20:end)),'ro');  hold on;      
        grid on;  axis equal;
        xlabel(['Re a_' num2str(iax)]); ylabel(['Im a_' num2str(iax)]);
        axis([-1 1 -1 1]);
      end
    case 3,  %  RPO stability properties
      load(name(ost(ist),:));
      [g, jac] = ksfmf5new(0, [a0; tend; ph], d, h, eye(N-2), [1 1]);
      disp(['T = ' num2str(tend) '  d = ' num2str(ph) ' |g| = ' num2str(norm(g))]);
      df = jac(1:N-2,1:N-2)+eye(N-2);
      [vdf, edf] = eig(df);  edf = diag(edf);
      [sedf, ie] = sort(abs(edf),1,'descend');
      edf = edf(ie);  vdf = vdf(:,ie);
      disp(sprintf('(%15.10f,%9.5f)\n',[abs(edf(1:10)) angle(edf(1:10))]'));
    end
  end

%% Plot FMs of RPOs in co-rotating frame
  clear;  d = 22;  N = 32;  x = d.*(-N:N)'./(2*N);  h = 0.25;
  load ks22uqo060a;
  [tt, aa] = ksfmedt(d, tend, a0, h, 1);  v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
  vv = [zeros(1,size(aa,2)); v; zeros(N+1,size(aa,2)); flipud(conj(v))];
  u = real(fft(vv)); u = [u; u(1,:)];
fig1 = figure('PaperPosition',[0.6345 6.345 20.3 15.23],...
    'PaperSize',[20.98 29.68],'Position',[5 400 1390 540]);
 ax1 = axes('Position',[0.04 0.09 0.13 0.85],'Parent',fig1);
  pcolor(x,tt,u'); shading flat;  caxis([-3 3]);  hold on; 
  title(sprintf('T = %6.2f  d = %9.6f',tend,ph));  
  xlabel('x'); ylabel('Time');

  ek = exp((-2i*pi/d).*ph./tend*(1:N/2-1)'*tt);  vr = ek.*v;
  e1 = exp(-2i*pi./tend*(1:N/2-1)'*tt);
  load ks22uss0a;  v2 = repmat(a0(1:2:end)+1i*a0(2:2:end),1,length(tt));
  load ks22uss3a;  v3 = repmat(a0(1:2:end)+1i*a0(2:2:end),1,length(tt));
  load ks22utw1a;  v4 = repmat(a0(1:2:end)+1i*a0(2:2:end),1,length(tt));
  v2r = e1.*v2;  v3r = e1.*v3;  v4r = e1.*v4;  v5r = -e1.*conj(v4);
 ax2 = axes('Position',[0.20 0.09 0.38 0.85],'Parent',fig1);
%  plot3(real(v(2,:)), imag(v(2,:)), real(v(3,:)), '.-');  hold on;
  plot3(real(vr(2,:)), imag(vr(2,:)), imag(vr(3,:)), '.-');
  grid on; hold on; 
  xlabel('Re u_2');  ylabel('Im u_2');  zlabel('Re u_3');
  plot3(real(v2r(2,:)), imag(v2r(2,:)), imag(v2r(3,:)), 'k.-');
  plot3(real(v3r(2,:)), imag(v3r(2,:)), real(v3r(3,:)), 'r.-'); 
%  plot3(real(v4r(2,:)), imag(v4r(2,:)), real(v4r(1,:)), 'g.-'); 
%  plot3(real(v5r(2,:)), imag(v5r(2,:)), real(v5r(1,:)), 'c.-'); 
  axis image;  view(-160,30);
 ax3 = axes('Position',[0.60 0.09 0.38 0.85],'Parent',fig1);
%  plot3(real(v(3,:)), imag(v(3,:)), real(v(2,:)), '.-');  hold on;
  plot3(real(vr(3,:)), imag(vr(3,:)), real(vr(2,:)), '.-');
  grid on; hold on;
  xlabel('Re u_3');  ylabel('Im u_3');  zlabel('Re u_2');
  plot3(real(v2r(3,:)), imag(v2r(3,:)), real(v2r(2,:)), 'k.-');
  plot3(real(v3r(3,:)), imag(v3r(3,:)), real(v3r(2,:)), 'r.-');
%  plot3(real(v4r(1,:)), imag(v4r(1,:)), real(v4r(2,:)), 'g.-');
%  plot3(real(v5r(1,:)), imag(v5r(1,:)), real(v5r(2,:)), 'c.-');
  axis image;  view(-160,30);
  
  
%% Bring RPOs to the standard reference point
  clear;  load ks22uqo110a;  d = 22;  N = 32;  x = d.*(1:N)'./N;  h = 0.1;
  ek = exp((-2i*pi/d).*ph.*(1:N/2-1)');
  [tt, aa] = ksfmedt(d, tend, a0, h, 1);
  ap = ek.*(aa(1:2:end,end)+1i*aa(2:2:end,end));
  norm(ap-(a0(1:2:end)+1i*a0(2:2:end))),
  ma = sqrt(aa(1:2:end,:).^2+aa(2:2:end,:).^2);
  figure(2); clf; plot(tt,ma(1:3,:),'.-'); grid on;
  [mm, im] = max(ma(1,:));  te = tt(im);
  hold on; plot(te, mm, 'ko');  a0 = aa(:,im);  
  pph = pi/2-angle(a0(1)+1i*a0(2));  
  ap = exp(1i*pph.*(1:N/2-1)').*(a0(1:2:end)+1i*a0(2:2:end));
  a0(1:2:end) = real(ap);  a0(2:2:end) = imag(ap);
  [tt, aa, da] = ksfmjedt(d, tend, a0, h);  eda = eig(da);
  ap = ek.*(aa(1:2:end)+1i*aa(2:2:end));
  norm(ap-(a0(1:2:end)+1i*a0(2:2:end))),
  cda = repmat(ek,1,N-2).*(da(1:2:end,:)+1i*da(2:2:end,:));
  da(1:2:end,:) = real(cda);  da(2:2:end,:) = imag(cda);  eda = eig(da);
  switch 2,
  case 1,
    [u,C,v] = svd(da-eye(N-2));  C = -v*u';  alpha = [1 1];%  C = eye(N-2);
    options = odeset('reltol',1e-7,'outputfcn',@ksdupoplot2);
    [s,y] = ode15s(@(t,a)ksfmf5new(t,a,d,h,C,alpha), [0 1000], ...
                   [a0; tend; ph], options);
  case 2,
    C = eye(N-2);  alpha = [1 1];
    options = optimset('Display','off','outputfcn',@ksoutfun,'MaxIter',2000, ...
              'MaxFunEvals',20000,'Jacobian','on','NonlEqnAlgorithm','lm','TolFun',1e-9);
    upo = fsolve(@(a)ksfmf5new(te,a,d,h,C,alpha),[a0;tend;ph],options);
  end, return;
  C = eye(N-2); alpha = [1 1]; eda = eig(C*da);
  [g, jac] = ksfmf5new(te, [a0;tend;ph], d, h, C, alpha);
  norm(g), % eja = eig(jac); [eda eja(1:end-2)],

%% Use time-and-space FT (after Lopez et al.)
  clear; load ks22seed140a;  d = 22;  N = 32;  x = d.*(1:N)'./N;  Nt = 64;  np = 4;
  h = tend/(Nt*np);  [tt, aa] = ksfmedt(d, tend, a0, h, np);
  figure(2); clf; plot(aa(1,:),aa(2,:),'.-');
  ek = exp((-2i*pi*ph/(d*Nt)).*((1:N/2-1)'*(0:Nt)));
  ap = ek.*(aa(1:2:end,:)+1i*aa(2:2:end,:));
  hold on; plot(real(ap(1,:)),imag(ap(1,:)),'r.-'); grid on;
  norm(ap(:,1)-ap(:,end)),
%  aa = zeros(N-2,Nt);  aa(1:2:end,:) = real(ap(:,1:end-1));  aa(2:2:end,:) = imag(ap(:,1:end-1));
%  for ii = 1:N/2,  figure(3); clf; plot(aa(2*ii-1,:),aa(2*ii,:),'.-'); grid on; pause; end,
  v = [zeros(1,Nt); ap(:,1:end-1); zeros(1,Nt); flipud(conj(ap(:,1:end-1)))];
  v = ifft(v,[],2);  u = real(fft2(v));  Nv = ifft2(u.^2);
  k = (2.*pi./d).*[0:N/2-1 0 -N/2+1:-1]';  L = k.^2 - k.^4;
  w = (2.*pi./tend).*[0:Nt/2-1 0 -Nt/2+1:-1];
  f = (repmat(L-1i*k*(ph/tend),1,Nt)+repmat(1i*w,N,1)).*v + 0.5i*repmat(k,1,Nt).*Nv;

%% Introduce intermediate points (ksfmms0 - multiple shooting)
  clear; load ks22uqo107a;  ni = [4 0];  d = 22;  N = 32;  h = 0.25;
  ni(2) = ceil(tend./(h.*(ni(1)+1)));
  [tt, aa] = ksfmedt(d, tend, a0, h, ni(2));
  ai = aa(:,1:end-1);  ai = [ai(:); tend-ni(1)*ni(2)*h; ph];
  hnew = 0.1;  ni(2) = round(ni(2)*h/hnew);  h = hnew;
  gi = ksfmfns(ai, d, h, ni);
%  load ks22seed123ai;
  options = optimset('Display','on', 'outputfcn',@ksoutfun,'MaxIter',2000, ...
      'MaxFunEvals',100000,'Jacobian','off','NonlEqnAlgorithm','lm','TolX',1e-12,'TolFun',1e-12);
  upo = fsolve(@(ai)ksfmfns(ai,d,h,ni),ai,options);

%% Evaluation of ksfmmint
  clear;  d = 22;  h = 0.25;  load ks22uqo040b;  tend = 50;
  [tt, aa] = ksfmmint(d, tend, a0, h, 1);
  g = aa-repmat(a0,1,size(aa,2));  gg = sum(g.^2)'.*h;  fa = aa;
  for ii = 1:size(aa,2),
    fa(:,ii) = ksfm(0, aa(:,ii), d); end
  fg = sum(fa.*g)';
  figure(1); clf; plot(tt,[gg fg],'.-'); grid on;

%% Analytic formula for minimum distance PSS
  clear;  syms a0 y0 y1 y0p y1p h dy dya s ys ysp a2 a3;
%  dy = (y1 - y0)./h;  dya = (y0 - a0)./h;
%  a2 = 3.*dy-y1p-2.*y0p;  a3 = y1p+y0p-2.*dy;
  ys = y0 + h*(y0p*s + a2*s^2 + a3*s^3);
  ysp = diff(ys,s);  cs = ysp*(ys-a0)/h^2;  collect(cs,s)

%% Explore various stabilising transformations for ksfmedt
  clear;  d = 22;  h = 0.25;  na = 30;  load ks22uqo042a;
  [tt, a, da] = ksfmjedt(d, tend, a0, h);
  k = (-2i*pi/d).*(1:na/2)';  ek = exp(k.*ph);
  c = ek.*(a(1:2:end,end)+1i*a(2:2:end,end));
  ap = a;  ap(1:2:end) = real(c);  ap(2:2:end) = imag(c);
  dc = repmat(ek,1,na).*(da(1:2:end,:)+1i*da(2:2:end,:));
  da(1:2:end,:) = real(dc);  da(2:2:end,:) = imag(dc);
  [U,T] = schur(da);  
  disp(sprintf(' %10.3f %10.3f %10.5f %10.5f %10.5f %10.5f\n',diag(T(1:6,1:6)),diag(T(2:7,1:6)))); 
  nr = 4;  r = U(:,nr);  la = T(1,1);
  figure(2); clf;
  for ir = (-4:0.1:4)./la,  ar = a0 + ir.*r;
    [tt, a] = ksfmjedt(d, tend, ar, h);
    c = ek.*(a(1:2:end,end)+1i*a(2:2:end,end));
    ap(1:2:end) = real(c);  ap(2:2:end) = imag(c);
    plot(ir,norm(ap-ar),'.'); hold on; end
  grid on;
%%  
  e = eig(da);  [ae, ie] = sort(abs(e),1,'descend');  disp(e(ie(1:6)));
  switch 2,
  case 1,  [R,T] = schur(da);  C = eye(na) - 2*R(:,2)*R(:,2)';
  case 2,  [U,S,V] = svd(da-eye(na));  C = -V*U';
  case 3,  [U,S,V] = svd(da-eye(na));  C = eye(na) - 2*U(:,1)*U(:,1)';
  end,
  es = eig(C*(da - eye(na)));  disp(es(1:6));
  
%% Stability of equilibria
  load ks22uss1a;  [h, dh] = ksfm(0, ksfmshift(a0), d);  [vdh, edh] = eig(dh);  edh = diag(edh);
  [sedh, ie] = sort(real(edh),1,'descend');  e1eig = edh(ie);  e1vec = vdh(:,ie);
  load ks22uss2a;  [h, dh] = ksfm(0, ksfmshift(a0), d);  [vdh, edh] = eig(dh);  edh = diag(edh);
  [sedh, ie] = sort(real(edh),1,'descend');  e2eig = edh(ie);  e2vec = vdh(:,ie);
  load ks22uss3b;  [h, dh] = ksfm(0, ksfmshift(a0), d);  [vdh, edh] = eig(dh);  edh = diag(edh);
  [sedh, ie] = sort(real(edh),1,'descend');  e3eig = edh(ie);  e3vec = vdh(:,ie);
  disp(' 1-wave equilibrium  2-wave equilibrium  3-wave equilibrium');
  disp([e1eig(1:10) e2eig(1:10) e3eig(1:10)]);
%   disp('1-wave equilibrium');
%   disp(sprintf(' (%9.5f,%9.5f)  (%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f, ... )',...
%     [real(e1eig(1)); imag(e1eig(1)); real(e1vec(1:6,1))]));
%   disp(sprintf(' (%9.5f,%9.5f)  (%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f, ... )',...
%     [real(e1eig(2)); imag(e1eig(2)); imag(e1vec(1:6,1))]));
%   disp(sprintf(' (%9.5f,%9.5f)  (%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f, ... )',...
%     [real(e1eig(3)); imag(e1eig(3)); real(e1vec(1:6,3))]));
%   disp(sprintf(' (%9.5f,%9.5f)  (%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f, ... )',...
%     [real(e1eig(4)); imag(e1eig(4)); imag(e1vec(1:6,3))]));
%   disp(sprintf(' (%9.5f,%9.5f)  (%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f, ... )',...
%     [real(e1eig(5)); imag(e1eig(5)); real(e1vec(1:6,5))]));
%   disp(sprintf(' (%9.5f,%9.5f)  (%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f, ... )',...
%     [real(e1eig(6)); imag(e1eig(6)); real(e1vec(1:6,6))]));
%   disp(sprintf(' (%9.5f,%9.5f)  (%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f, ... )',...
%     [real(e1eig(7)); imag(e1eig(7)); imag(e1vec(1:6,6))]));
%   disp(sprintf(' (%9.5f,%9.5f)  (%7.3f,%7.3f,%7.3f,%7.3f,%7.3f,%7.3f, ... )\n',...
%     [real(e1eig(8:10)'); imag(e1eig(8:10)'); real(e1vec(1:6,8:10))]));
  
%% Stability of equilibria
  load ks22uss1a;  [h, dh] = ksfm(0, a0, d);  [vdh, edh] = eig(dh);  edh = diag(edh);
  [sedh, ie] = sort(real(edh),1,'descend');  e1eig = edh(ie);  e1vec = vdh(:,ie);
  load ks22uss2a;  [h, dh] = ksfm(0, a0, d);  [vdh, edh] = eig(dh);  edh = diag(edh);
  [sedh, ie] = sort(real(edh),1,'descend');  e2eig = edh(ie);  e2vec = vdh(:,ie);
  load ks22uss3b;  [h, dh] = ksfm(0, a0, d);  [vdh, edh] = eig(dh);  edh = diag(edh);
  [sedh, ie] = sort(real(edh),1,'descend');  e3eig = edh(ie);  e3vec = vdh(:,ie);
  disp(' 1-wave equilibrium  2-wave equilibrium  3-wave equilibrium');
  disp([e1eig(1:10) e2eig(1:10) e3eig(1:10)]);

%% Evolution of orbits starting in the unstable manifolds of the equilibria
  load ks22uss1a;  tend = 150;  np = 3;
  a1 = a0 + 1e-3.*real(e1vec(:,1)); 
  [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
  figure(1); set(gcf,'pos',[100 500 800 400]); clf;
  ax1 = axes('pos',[0.06 0.10 0.20 0.80]); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat; title('1-wave, plane 1');

  for angl = (0:.05:2)*pi,
    a1 = a0 + 1e-3.*(cos(angl).*real(e1vec(:,1))+sin(angl).*imag(e1vec(:,1)));
    [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
    pcolor(x,tt,uu'); caxis([-3 3]); shading flat; title(['\alpha = ' num2str(angl)]);
    pause; end,

  a1 = a0 + 1e-3.*real(e1vec(:,3));
  [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
  ax2 = axes('pos',[0.30 0.10 0.20 0.80]); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat; title('1-wave, plane 2');

  for angl = (0:.05:2)*pi,
    a1 = a0 + 1e-3.*(cos(angl).*real(e1vec(:,3))+sin(angl).*imag(e1vec(:,3)));
    [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
    pcolor(x,tt,uu'); caxis([-3 3]); shading flat; title(['\alpha = ' num2str(angl)]);
    pause; end,

  load ks22uss2a;  tend = 150;  np = 2;
  a1 = a0 + 1e-3.*imag(e2vec(:,1)); 
  [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
  ax3 = axes('pos',[0.54 0.10 0.20 0.80]); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat; title('2-wave equilibrium');

  for angl = (0:.05:2)*pi,
    a1 = a0 + 1e-3.*(cos(angl).*real(e2vec(:,1))+sin(angl).*imag(e2vec(:,1)));
    [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
    pcolor(x,tt,uu'); caxis([-3 3]); shading flat; title(['\alpha = ' num2str(angl)]);
    pause; end,

  load ks22uss3b;  tend = 150;  np = 2;
  a1 = a0 + 1e-3.*e3vec(:,1);
  [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
  ax4 = axes('pos',[0.78 0.10 0.20 0.80]); pcolor(x,tt,uu'); caxis([-3 3]);
    shading flat; title('3-wave equilibrium');
    
  for angl = (0:.05:2)*pi,
    a1 = a0 + 1e-3.*(cos(angl).*e3vec(:,1)+sin(angl).*e3vec(:,2));
    [tt, aa] = ksfmedt(d, tend, a1, h, np);  [x, uu] = ksfm2real(aa, d, 64);
    pcolor(x,tt,uu'); caxis([-3 3]); shading flat; title(['\alpha = ' num2str(angl)]);
    pause; end,

%% Improve accuracy of equilibria and save results in kse22orbits
% eq(k).a; eq(k).eig; eq(k).evec;
% save kse22orbits eq tw rpo L h
  clear;  load kse22orbits;
  L = 22;  N = 32;  % number of modes (number of ODEs = 2*N-2)
  names = {'ks22uss1a'; 'ks22uss2a'; 'ks22uss3a'};
  for k = 1:length(names),
    load(names{k});
    if length(a0) < 2*N-2, n = 2*N-2-length(a0); end
    ae = [ksfmshift(a0); zeros(n,1)];
    options = optimset('Display','iter','Jacobian','on','TolX',1e-12,'TolFun',1e-15);
    ae = fsolve(@(a)ksfm(0,a,L),ae,options);  ae = ksfmshift(ae);
    fae = ksfm(0,ae,L);  disp(['|f(ae)| = ' num2str(norm(fae))]);
%    [tt,aa] = ksfmedt(L, 500, ae, 0.1, 10); [x, uu] = ksfm2real(aa, L);
%    figure(1); clf; pcolor(x,tt,uu'); tall(1.1); shading flat; caxis([-3 3]); pause; 
    [hh, dh] = ksfm(0, ae, L);  [vdh, edh] = eig(dh);  edh = diag(edh);
    [sedh, ie] = sort(real(edh),1,'descend');
    eq(k).a = ae;  eq(k).eig = edh(ie);  eq(k).evec = vdh(:,ie);
  end
  save kse22orbits eq L

%% Improve accuracy of relative equilibria and save results in kse22orbits
% re(k).a;  re(k).c;  re(k).eig;  re(k).evec;
  clear;  L = 22;  N = 32;  % number of modes (number of ODEs = 2*N-2)
  names = {'ks22utw1a'; 'ks22utw2a'};
  for k = 1:length(names),
    load(names{k});
    if length(a0) < 2*N-2, n = 2*N-2-length(a0);  
      are = [ksfmshift(a0); zeros(n,1)];  
    else are = ksfmshift(a0); end
    if c < 0; are(1:2:end) = -are(1:2:end); c = -c; end
    options = optimset('Display','iter','Jacobian','on','TolX',1e-14,'TolFun',1e-15,'DerivativeCheck','on');
    aa = fsolve(@(a)ksfmtr(a,L),[are; c],options);
    re(k).a = ksfmshift(aa(1:end-1));  re(k).c = aa(end);
    [f,df] = ksfmtr([re(k).a;re(k).c],L); disp(['|f(re)| = ' num2str(norm(f))]);
%   [tt,aa] = ksfmedt(L, 100, re(1).a, 0.1, 2); [x, uu] = ksfm2real(aa, L);
%   figure(1); clf; pcolor(x,tt,uu'); tall(1.1); shading flat; caxis([-3 3]); pause; 
    [vdf, edf] = eig(df(1:end-1,1:end-1));  edf = diag(edf);
    [sedf, ie] = sort(real(edf),1,'descend');
    re(k).eig = edf(ie);  re(k).evec = vdf(:,ie); end
  save kse22orbits re -append

%% Improve accuracy of RPOs and save results in kse22orbits
% rpo(k).a;  rpo(k).T;  rpo(k).d;  rpo(k).eig;  rpo(k).evec;  
  clear;  load kse22orbits;
  L = 22;  N = 32;  % number of modes (number of ODEs = 2*N-2)
  h = 0.1;          % integration stepsize for ksfmedt
  switch 2,
  case 1, load ks22names;  iname = 102; load(name(iname,:),'a0','tend','ph');
  case 2, load('ks22uqo183b','a0','tend','ph');
  case 3, ipo = 99;  a0 = rpo(ipo).a;  tend = rpo(ipo).T;  ph = rpo(ipo).d; end
  disp(['Period = ' num2str(tend) '  Phase = ' num2str(ph)]);
  if 1,  %%% use multiple shooting  
    ni = [6 0];  ho = 0.25;  ni(2) = 2*ceil(tend./(ho.*(ni(1)+1))/2);
    disp(['Using multiple shooting: ' num2str(ni)]);
    [tt, aa] = ksfmedt(L, tend, a0, ho, ni(2));
    if (size(aa,1) < N-2), aa = [aa; zeros(N-2-size(aa,1),size(aa,2))]; end
    ai = aa(:,1:end-1);  ai = [ai(:); tend-ni(1)*ni(2)*ho; ph];
    ni(2) = round(ni(2)*ho/h);  gi = ksfmfns(ai, L, h, ni);
    disp(['|gms| = ' num2str(norm(gi))]);
    options = optimset('Display','off', 'outputfcn',@ksoutfun2,'MaxIter',2000, ...
      'MaxFunEvals',100000,'Jacobian','off','NonlEqnAlgorithm','lm','TolX',1e-12,'TolFun',1e-12);
    upo = fsolve(@(ai)ksfmfns(ai,L,h,ni),ai,options);
    a0 = upo(1:(length(upo)-2)/(ni(1)+1));  
    tend = upo(end-1)+ni(1)*ni(2)*h;  ph = upo(end);
    disp(['Period = ' num2str(tend) '  Phase = ' num2str(ph)]);  end
  if length(a0) < 2*N-2, n = 2*N-2-length(a0); a0 = [ksfmshift(a0); zeros(n,1)];
  else a0 = ksfmshift(a0); end
  if ph > L/2, ph = ph-L; elseif ph < -L/2, ph = ph+L; end
  if ph < 0, ph = -ph; a0(1:2:end) = -a0(1:2:end); end
  options = optimset('Display','off','outputfcn',@ksoutfun2,'MaxIter',5000, ...
    'MaxFunEvals',20000,'Jacobian','on','NonlEqnAlgorithm','lm','TolX',1e-12,'TolFun',1e-12);
  [upo,fval,eflag] = fsolve(@(a)ksfmf5(0,a,L,h,eye(2*N),[1 1]),[a0;tend;ph],options);
  a0 = upo(1:end-2); tend = upo(end-1); ph = upo(end);
  for ic = 1:20, format short g;
    [g, jac] = ksfmf5(0,[a0;tend;ph],L,h,eye(2*N),[1 1]);
    disp(['|g|=' num2str(norm(g(1:end-2)))]);  de = eig(jac(1:end-2,1:end-2)+eye(2*N-2)); 
    [se,ie]=sort(abs(de),'descend');  disp(de(ie(1:10)));
    [u,C,v] = svd(jac);  C = -v*u';  alpha = [1e0 1e0];
%    [u,C,v] = svd(jac(1:end-2,1:end-2)); C = blkdiag(-v*u',eye(2)); format short;
    options = odeset('abstol',1e-15,'reltol',2.3e-14,'outputfcn',@ksdupoplot4);
    [s,y] = ode15s(@(t,a)ksfmf5(t,a,L,h,C,alpha), [0 100000], [a0; tend; ph], options);
    a0 = y(end-1,1:end-2)'; tend = y(end-1,end-1), ph = y(end-1,end),
    but = input('Stop? ','s');  if but == 'y', break; end, end
% Determine if this orbit already exists in rpo
  if exist('rpo'), ipo = max(find([rpo.T] < tend+1e-4)); add = 1;
   if isempty(ipo), ipo = 0; add = 1;
   elseif abs(rpo(ipo).T-tend) < 1e-4,  % exists
     disp(['This RPO appears similar to the existing one: ipo = ' num2str(ipo)]);
     disp(num2str([rpo(ipo).T tend; rpo(ipo).d ph])); add = 0; end
  else ipo = 0;  add = 1;  rpo = []; end
  if add == 1, ipo = ipo+1;
    for iit = size(rpo,2):-1:ipo,  rpo(iit+1) = rpo(iit); end
    rpo(ipo).a = a0;  rpo(ipo).T = tend;  rpo(ipo).d = ph;
    [f,df] = ksfmf5(0,[rpo(ipo).a;rpo(ipo).T;rpo(ipo).d],L,h,eye(2*N),[1 1]); 
    disp(['|f|=' num2str(norm(f(1:end-2)))]);
    [vdf,edf] = eig(df(1:end-2,1:end-2)+eye(2*N-2)); edf = diag(edf);
    [sedf, ie] = sort(abs(edf),1,'descend');
    rpo(ipo).eig = edf(ie);  rpo(ipo).evec = vdf(:,ie);
    disp(['Added orbit ' num2str(ipo) ' of ' num2str(size(rpo,2))]);
    disp(['Period = ' num2str(tend) '  Phase = ' num2str(ph)]); 
    save kse22orbits rpo h -append
  end
  
%% Explore unstable manifolds of equilibria (3-wave)
  clear;  load kse22orbits;  k = 3;  h = 0.1;  delta = 1e-4;  tend = 200;
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) eq(k).evec(:,4)]);
  figure(1); clf;  figure(2); clf; 
  for phi = 0.2532017200:0.0000000001:0.2532017208, %(0:0.1:6)*pi./3,
    a0 = eq(k).a + delta.*(v(:,1).*cos(phi)+v(:,2).*sin(phi));
    [tt, aa] = ksfmedt(L, tend, a0, h, 2);
    [x, uu] = ksfm2real(aa, L, 64);  av = v'*aa;
    figure(1);  pcolor(x,tt,uu'); caxis([-3 3]); shading flat;  hold on;
    title(['\phi = ' sprintf('%13.10f',phi)]);
    figure(2);  plot3(av(1,:),av(2,:),av(3,:),'-'); 
    hold on; grid on; axis equal; view(0,90); pause;
  end,

%% Explore unstable manifolds of equilibria (2-wave)
  clear;  load kse22orbits;  k = 2;  h = 0.1;  tend = 400;
  ere = real(eq(k).eig(1));  eim = imag(eq(k).eig(1));
  period = 2*pi/eim;  scale = exp(ere*period);
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
  figure(1); clf;  figure(2); clf; 
  for delta = 0.2535144:.0000001:0.2535145, %0:0.1:ere*period,
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmedt(L, tend, a0, h, 5);
    [x, uu] = ksfm2real(aa, L, 64);  av = v'*aa;
    figure(1);  pcolor(x,tt,uu'); caxis([-3 3]); shading flat;  hold on;
    title(['\delta = ' sprintf('%12.9f',delta)]);
    figure(2);
    plot3(av(1,:),av(2,:),av(3,:),'-'); hold on; grid on; axis equal; pause;
 %   plot(av(1,:),av(2,:),'-'); hold on; grid on; axis equal; pause;
 %   plot3(tt,av(1,:),av(2,:),'.-'); hold on; grid on; view(0,0); pause;
  end,

%% Explore unstable manifolds of equilibria (1-wave, plane 1,2: p = 1,3)
  clear;  load kse22orbits;  k = 1;  p = 3;  h = 0.1;  tend = 200;
  ere = real(eq(k).eig(p));  eim = imag(eq(k).eig(p));
  period = 2*pi/eim;  scale = exp(ere*period);
  v = gsorth([real(eq(k).evec(:,p)) imag(eq(k).evec(:,p)) real(eq(k).evec(:,6))]);
  figure(1); clf;  figure(2); clf; 
  for delta = 0:0.05:ere*period,
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmedt(L, tend-delta/ere, a0, h, 2);
    [x, uu] = ksfm2real(aa, L, 64);  av = v'*aa;
    figure(1);  pcolor(x,tt,uu'); caxis([-3 3]); shading flat;
    title(['\delta = ' num2str(delta)]);
    figure(2);
    plot3(av(1,:),av(2,:),av(3,:),'-'); hold on; grid on; axis equal; view(-30,30); pause;
 %   plot(av(1,:),av(2,:),'-'); hold on; grid on; axis equal; pause;
  end,

%% Plot equilibria and travelling waves in (u,u_x,u_xx) coordinates
  clear;  load kse22orbits; colr = ['b';'r';'k'];
  figure(1); set(gcf,'pos',[420 500 560 450]); clf;
  ax1 = axes('Position',[0.05 0.0828 0.788 0.907]);
  axis(ax1,[-2.597 2.597 -3.002 1.459 -2.804 2.804]);
  xlabel(ax1,'u'); ylabel(ax1,'u_x'); zlabel(ax1,'u_{xx}','rotat',0);
  view(ax1,[25 25]); grid(ax1,'on'); hold(ax1,'all');
  for k = 1:3, [x,u,ux,uxx] = ksfm2real(eq(k).a,L,128);
    plot3(u,ux,uxx,['.-' colr(k)]); hold on; end
  grid on; axis equal; set(gca,'pos',[0.03 0.05 0.9 0.94]); view(25,25);
  lg1 = legend(ax1,{'1-wave','2-wave','3-wave'});
  set(lg1,'pos',[0.75 0.3 0.18 0.14]);  set(ax1,'pos',[0.02 0.0828 0.788 0.907]);
  set(get(gca,'xlabel'),'pos',[17 -39 15.2]);  set(get(gca,'ylabel'),'pos',[20 -37 15.5]);
  set(get(gca,'zlabel'),'pos',[-1.9 -5.8 1.6]);
  
%% PCA for ksfm
  clear; 
  d = 22;  h = 0.25;  N = 20;  x = d.*(1:N)'./N;  np = 10;  tsampl = 2020;
  a0 = zeros(N-2,1); randn('seed',22000001);  a0(1:6) = 0.2*randn(6,1);
  [tt, aa, a0p, ap] = ksfmedt(d, tsampl, a0, h, np);  na = size(aa,2);
  v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
  v = [zeros(1,size(aa,2)); v; zeros(1,size(aa,2)); flipud(conj(v))];
  u = real(fft(v));
  figure(1); set(gcf,'pos',[5 35 250 910]); clf; 
  pcolor(x,tt,u'); shading flat; tall(1.1); hold on;

%%
  X = aa(:,6:end);  % sample set;
  mX = mean(X')';  B = X - repmat(mX,1,size(X,2));
  C = (B*B')./size(X,2);  [V,D] = eig(C);  d = diag(D);
  [sd, id] = sort(d,1,'descend');  V = V(:,id);  D = diag(sd);
  cumE = cumsum(sd)./sum(sd);  nP = min(find(cumE > 0.99));
  Vn = V(:,1:nP);  v = Vn(1:2:end,:) + 1i*Vn(2:2:end,:);
  v = [zeros(1,size(Vn,2)); v; zeros(1,size(Vn,2)); flipud(conj(v))];
  u = real(fft(v));
  figure(2); clf;  plot(x,u,'.-');

%% Plot RPOs in real space
  clear;  load kse22orbits; nax = 6;  np = 5;  te = 200;
% for ist = 1:nst,  ost = 1:nst;
% for ist = 1:nax, ost = 2:6, % ost = [3 22 50 57 95];
  for ist = 1:nax,  ost = [2 23 52 58 113 151]; %ost = [73:78];
    ifig = floor(ist/nax)+1;  iorb = mod(ist-1,nax)+1;
    if iorb == 1,  figure(ifig); clf; orient landscape; 
%      set(gcf,'pos',[35 200 1100 740]);
      set(gcf,'pos',[35 200 1000 640]);
      hax = subplots(1,nax,[0.03 0.03 0.07 0.1],[0.015 0 0 0]); end

    a0 = rpo(ost(ist)).a;  tend = rpo(ost(ist)).T;  ph = rpo(ost(ist)).d;
    N = length(a0)+2;  ek = exp((2i*pi/L).*ph.*(1:N/2-1)');
    [tt, aa] = ksfmedt(L, tend, a0, h, 1);
    tti = tt(1:end-1);  aai = aa(:,1:end-1);  ne = ceil(te./tend);  aav = aa;
    for ie = 1:ne-1,
      vi = (aav(1:2:end,:)+1i*aav(2:2:end,:)).*repmat(ek,1,size(aav,2));
      aav(1:2:end,:) = real(vi);  aav(2:2:end,:) = imag(vi);
      aai = [aai aav(:,1:end-1)];  tti = [tti tt(1:end-1)+ie*tend]; end
    [x, uui] = ksfm2real(aai, L);
          
    axes(hax(iorb)); pcolor(x,tti,uui'); shading flat; caxis([-3 3]);  hold on;
    plot(x([1 end])*ones(1,ne),[tend;tend]*(1:ne),'w-');
    plot(mod([ph;ph]*(1:ne-1),L)-L/2,[(1:ne-1);(2:ne)]*tend,'w-');
    title([sprintf('%3d: T=%6.2f ',ost(ist),tend) '\Delta=' sprintf('%7.3f',ph)]);
    set(gca,'ticklength',0.7*get(gca,'ticklength'),'tickdir','out',...
        'ylim',[0 te],'ytick',20:20:te,'box','off');  xlabel('x');
    if iorb > 1, set(gca,'yticklabel',[]); end
    if iorb == nax,
      hax(nax+1) = axes('pos',[.04 .9 .9 .07]);
      tit1 = 'Kuramoto-Sivashinsky: u_t = -uu_x - u_{xx} - u_{xxxx};  ';
      tit2 = ['x \in [-L/2, L/2];    BC: u(x+L,t) = u(x,t);' sprintf('   L = %4.1f;    ',L)];
      tit3 = 'Solutions of the form: u(x+\Delta,T) = u(x,0)';
      text(0.1,.9,[tit1 tit2 tit3]);  set(hax(nax+1),'visi','off');
      axes(hax(1)); ylabel('t','rotat',0); end
  end
  
%% Plot individual RPOs with the first 4 FMs in polar coordinates
  clear;  load kse22orbits;  te = 200;  np = 1;  refl = 0;
  for ist = 1:156, ost = 1:156; %ost = [1 2 5 10 50]; %ost = 20:5:40;
    a0 = rpo(ost(ist)).a;  tend = rpo(ost(ist)).T;  ph = rpo(ost(ist)).d;
    if refl, a0r = a0;  a0r(1:2:end) = -a0r(1:2:end);  phr = -ph; end
    eig = rpo(ost(ist)).eig(1:8);
    N = length(a0)+2;  ek = exp((2i*pi/L).*ph.*(1:N/2-1)');
    [tt, aa] = ksfmedt(L, tend, a0, h, 1);
    if refl, [tr, ar] = ksfmedt(L, tend, a0r, h, 1); end,
    if 1, % shift orbit in time
      [maa, ima] = max(sum(aa.^2)');  %disp([ima size(aa,2)]);
      if ima > 1, tt = [tt(ima:end) tt(2:ima)+tend]-tt(ima);
        aas = aa(:,2:ima);
        vi = (aas(1:2:end,:)+1i*aas(2:2:end,:)).*repmat(ek,1,ima-1);
        aas(1:2:end,:) = real(vi);  aas(2:2:end,:) = imag(vi);
        aa = [aa(:,ima:end) aas]; end, end
    if 0, phase = 5.351, % shift orbit's phase
      vi = (aa(1:2:end,:) + 1i*aa(2:2:end,:)).*...
           repmat(exp((2i*pi/L).*phase.*(1:N/2-1)'),1,length(tt));
      aa(1:2:end,:) = real(vi);  aa(2:2:end,:) = imag(vi);
    end
    tti = tt(1:end-1);  aai = aa(:,1:end-1);  ne = ceil(te./tend);  aav = aa;
    for ie = 1:ne-1,
      vi = (aav(1:2:end,:)+1i*aav(2:2:end,:)).*repmat(ek,1,size(aav,2));
      aav(1:2:end,:) = real(vi);  aav(2:2:end,:) = imag(vi);
      aai = [aai aav(:,1:end-1)];  tti = [tti tt(1:end-1)+ie*tend]; end
    [x, uui] = ksfm2real(aai, L);
    
    fig1 = figure('PaperPosition',[1.0 1.0 26 18.5],'PaperOrient','port',...
         'PaperSize',[20.98 29.68],'Position',[35 200 1000 740]);
    ax(1) = axes('Position',[0.06 0.08 0.15 0.86],'Parent',fig1);
    ax(2) = axes('Position',[0.24 0.52 0.36 0.42],'Parent',fig1);
    ax(3) = axes('Position',[0.64 0.52 0.36 0.42],'Parent',fig1);
    ax(4) = axes('Position',[0.24 0.03 0.36 0.42],'Parent',fig1);
    ax(5) = axes('Position',[0.64 0.03 0.36 0.42],'Parent',fig1);
    ax(6) = axes('Position',[0.22 0.08 0.02 0.86],'Parent',fig1);
    ax(7) = axes('Position',[0.60 0.56 0.01 0.02],'Parent',fig1);
    axes(ax(1));  cmap = jet;
    pcolor(x,tti,uui'); caxis([-3 3]); shading flat; hold on;
    plot(x([1 end])*ones(1,ne),[tend;tend]*(1:ne),'w-');
    plot(mod([ph;ph]*(1:ne-1),L)-L/2,[(1:ne-1);(2:ne)]*tend,'w-');
    title([sprintf('%3d: T=%6.2f ',ost(ist),tend) '\Delta=' sprintf('%7.3f',ph)]); 
    set(gca,'ticklength',0.7*get(gca,'ticklength'),'tickdir','out',...
      'ylim',[0 te],'ytick',20:20:te,'box','off');  xlabel('x'); ylabel('Time');

    for iax = 1:4,   axes(ax(iax+1));
      vi = aa(2*iax-1,:) + 1i*aa(2*iax,:);  vim = max(abs(vi));
      polar(0,vim,'w.'); hold on;
      polarcol(angle(vi),abs(vi),tt,cmap);
      if refl, hold on; vr = ar(2*iax-1,:) + 1i*ar(2*iax,:);
        polarcol(angle(vr),abs(vr),tr,cmap); end
      title(sprintf('Mode %d',iax),'fontsize',14); end
    axes(ax(6)); plotcol(tt.*0,tt,tt,cmap); set(gca,'vis','off','ylim',[0 te]);
    axes(ax(7));
    text(1,0,sprintf('(%8.5f,%8.3f)\n',[log(abs(eig))./tend angle(eig)./pi]'),'vert','top','horiz','cen');
    text(0,1,'$(\log~|\lambda|/T,~\arg~(\lambda)/\pi)$','interp','latex','horiz','cen');
    set(gca,'vis','off');
    print(gcf,'-dpng',sprintf('ks22rpo%06.2f_%06.3f.png',tend,ph)); delete(fig1);
  end
  
%% Determining proximity of chaotic orbit to RPOs
  clear;  load kse22orbits;  h = 0.1;  N = 32;  nph = 120;
  if 0,  np = 2;  tpre = 100;  tsampl = 1000;
    a0 = zeros(N-2,1); randn('seed',22000025);  a0(1:10) = 0.2*randn(10,1);
    [tt, aa] = ksfmedt(L, tpre, a0, h, 0);  a0 = aa;
    [tt, aa] = ksfmedt(L, tsampl, a0, h, np);
  else, load ks22prox04; end
  figure(1); clf; subplot(2,1,1); [x, uu] = ksfm2real(aa, L, 64);
  pcolor(tt,x,uu); caxis([-3 3]); shading flat; wide(1.2); pause;
  subplot(2,1,2); wide(1.2);
  va = aa(1:2:end,:) + 1i*aa(2:2:end,:);  na = size(aa,2);
  for ieq = 1:3,                      % Proximity to Equilibria
    arpo = eq(ieq).a;  vr = arpo(1:2:N-2) + 1i*arpo(2:2:N-2);
    for ib = 1:na,  phase = exp((-2i*pi/nph).*(1:N/2-1)'*(1:nph));
      dv = repmat(vr,1,nph).*phase - repmat(va(:,ib),1,nph);
      dd = sqrt(sum(real(dv).^2+imag(dv).^2));
      [mdd, idd] = min(dd);  prox(ieq,ib) = mdd;
      tprox(ieq,ib) = 0;  phprox(ieq,ib) = (idd/nph-0.5)*L;  end;
    plot(tt,prox(ieq,:),'.-'); hold all;
    disp([ieq]); pause; end
  for itw = 1:2, for irfl = -1:0,     % Proximity to Traveling Waves
    arpo = re(itw).a; if irfl == 0, arpo(1:2:end) = -arpo(1:2:end); end
    vr = arpo(1:2:N-2) + 1i*arpo(2:2:N-2);
    for ib = 1:na,  phase = exp((-2i*pi/nph).*(1:N/2-1)'*(1:nph));
      dv = repmat(vr,1,nph).*phase - repmat(va(:,ib),1,nph);
      dd = sqrt(sum(real(dv).^2+imag(dv).^2));  [mdd, idd] = min(dd);
      prox(2*itw+irfl+3,ib) = mdd;
      tprox(2*itw+irfl+3,ib) = 0;  
      phprox(2*itw+irfl+3,ib) = (idd/nph-0.5)*L;  end;
    plot(tt,prox(2*itw+irfl+3,:),'.-'); hold all;
    disp([itw irfl]); pause; end, end  
  for irpo = 61:70, for irfl = -1:0, tic;% Proximity to RPOs
    iprox = 2*irpo+irfl+7;  arpo = rpo(irpo).a;  
    proxTrpo(iprox) = rpo(irpo).T;
    proxDrpo(iprox) = rpo(irpo).d;
    if irfl == 0, arpo(1:2:end) = -arpo(1:2:end); 
      proxDrpo(iprox) = -proxDrpo(iprox); end
    [tr, ar] = ksfmedt(L, rpo(irpo).T, arpo, h, 1);  nr = size(ar,2);
    vr = ar(1:2:N-2,:) + 1i*ar(2:2:N-2,:);  dd = zeros(nph,nr);
    if irfl == -1, phase = exp((-2i*pi*rpo(irpo).d./L).*(1:N/2-1)');
    else, phase = exp((2i*pi*rpo(irpo).d./L).*(1:N/2-1)'); end
    disp(sqrt(sum(abs(vr(:,end).*phase - vr(:,1)).^2)));
    for ib = 1:na,
      for iph = 1:nph,  phase = exp((-2i*pi/nph*iph).*(1:N/2-1)');
        dv = vr.*repmat(phase,1,nr) - repmat(va(:,ib),1,nr);
        dd(iph,:) = sqrt(sum(real(dv).^2+imag(dv).^2)); end
      [mdd, idd] = min(dd(:));  [idr,idc] = ind2sub(size(dd),idd);
      prox(iprox,ib) = dd(idr,idc);  tprox(iprox,ib) = tr(idc);     
      phprox(iprox,ib) = idr/nph*L;  end;
    plot(tt,prox(iprox,:),'.-'); hold all;
    disp([irpo irfl iprox proxTrpo(iprox) proxDrpo(iprox)]); 
    drawnow; toc; end, end

%% Determining proximity of RPOs to equilibria and each other
  clear;  load kse22orbits;  h = 0.1;  N = 32;  nph = 120;
  for ir1 = 2:2,
    [tr1, ar1] = ksfmedt(L, rpo(ir1).T, rpo(ir1).a, h, 1);
    figure(1); clf; subplot(2,1,1); [x, uu] = ksfm2real(ar1, L, 64);
    pcolor(tr1,x,uu); caxis([-3 3]); shading flat; wide(1.2);
    subplot(2,1,2); wide(1.2);
    vr1 = ar1(1:2:N-2,:) + 1i*ar1(2:2:N-2,:);  nr1 = size(ar1,2);
    disp(sqrt(sum(abs(vr1(:,end).*exp((-2i*pi*rpo(ir1).d./L).*(1:N/2-1)') - vr1(:,1)).^2)));
    pause;
    for ieq = 1:3, prox = [];             % Proximity to Equilibria
      aeq = eq(ieq).a;  veq = aeq(1:2:N-2) + 1i*aeq(2:2:N-2);
      for ib = 1:nr1,  phase = exp((-2i*pi/nph).*(1:N/2-1)'*(1:nph));
        dv = repmat(veq,1,nph).*phase - repmat(vr1(:,ib),1,nph);
        dd = sqrt(sum(real(dv).^2+imag(dv).^2)); [mdd, idd] = min(dd);
        prox(ib) = mdd;  phprox(ib) = (idd/nph-0.5)*L; end
      [mprox,iprox] = min(prox);
      plot(tr1,prox,'.-'); hold all; set(gca,'xlim',[0 rpo(ir1).T]);
      plot(tr1(iprox),mprox,'ko');
      disp([ieq mprox tr1(iprox) phprox(iprox)]); pause; end
    for itw = 1:2, for irfl = -1:0,     % Proximity to Traveling Waves
      atw = re(itw).a; if irfl == 0, atw(1:2:end) = -atw(1:2:end); end
      vtw = atw(1:2:N-2) + 1i*atw(2:2:N-2);  prox = [];
      for ib = 1:nr1,  phase = exp((-2i*pi/nph).*(1:N/2-1)'*(1:nph));
        dv = repmat(vtw,1,nph).*phase - repmat(vr1(:,ib),1,nph);
        dd = sqrt(sum(real(dv).^2+imag(dv).^2));  [mdd, idd] = min(dd);
      prox(ib) = mdd;  phprox(ib) = (idd/nph-0.5)*L;  end;
      [mprox,iprox] = min(prox);      
      plot(tr1,prox,'.-'); hold all; set(gca,'xlim',[0 rpo(ir1).T]);
      plot(tr1(iprox),mprox,'ko');
      disp([itw irfl mprox tr1(iprox) phprox(iprox)]); pause; end, end  
  end
%%  
  for irpo = 33:33, for irfl = -1:0, tic;% Proximity to RPOs
    iprox = 2*irpo+irfl+7;  arpo = rpo(irpo).a;  
    proxTrpo(iprox) = rpo(irpo).T;
    proxDrpo(iprox) = rpo(irpo).d;
    if irfl == 0, arpo(1:2:end) = -arpo(1:2:end); 
      proxDrpo(iprox) = -proxDrpo(iprox); end
    [tr, ar] = ksfmedt(L, rpo(irpo).T, arpo, h, 1);  nr = size(ar,2);
    vr = ar(1:2:N-2,:) + 1i*ar(2:2:N-2,:);  dd = zeros(nph,nr);
    for ib = 1:na,
      for iph = 1:nph,  phase = exp((-2i*pi/nph*iph).*(1:N/2-1)');
        dv = vr.*repmat(phase,1,nr) - repmat(va(:,ib),1,nr);
        dd(iph,:) = sqrt(sum(real(dv).^2+imag(dv).^2)); end
      [mdd, idd] = min(dd(:));  [idr,idc] = ind2sub(size(dd),idd);
      prox(iprox,ib) = dd(idr,idc);  tprox(iprox,ib) = tr(idc);     
      phprox(iprox,ib) = idr/nph*L;  end;
    plot(tt,prox(iprox,:),'.-'); hold all;
    disp([irpo irfl iprox proxTrpo(iprox) proxDrpo(iprox)]); 
    drawnow; toc; end, end

%% Plotting E2 -> E2 orbits with the unstable manifold of E2
  clear;  load kse22orbits;  irpo = [15 17 19 29 35 45 48 63]; av = [];
  tend = 150; rcol = ['b';'r';'k';'g';'c';'m';'y'];  figure(1); clf;
  ere = real(eq(2).eig(1));  period = 2*pi/imag(eq(2).eig(1));
  v = gsorth([real(eq(2).evec(:,1)) imag(eq(2).evec(:,1)) real(eq(2).evec(:,7))]);
  for delta = (0:0.02:1.02).*ere*period,
    a0 = eq(2).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmedt(L, tend, a0, h, 2); av = [av; v'*aa]; end,
  plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','k-');
  hold on; grid on; axis equal; view(10,20);
  for ir = 1:8,
    [tr, ar] = ksfmedt(L, rpo(irpo(ir)).T, rpo(irpo(ir)).a, h, 1);
    % Translate RPO in time and phase to the point closest to E2
    vr = ar(1:2:end,:) + 1i*ar(2:2:end,:);
    veq = eq(2).a(1:2:end) + 1i*eq(2).a(2:2:end);
    vrs = vr.*exp(1i*(1:31)'*(pi/4-angle(vr(2,:))/2));
    [mv, itv] = min(abs(vrs(2,:)-veq(2)));
    [tr, ar] = shiftrpo(tr, ar, L, rpo(irpo(ir)).T, rpo(irpo(ir)).d, itv);
    vr = ar(1:2:end,:) + 1i*ar(2:2:end,:);  phase = pi/4 - angle(vr(2,1))/2;
    vr = vr.*repmat(exp(1i.*phase.*(1:size(vr,1))'),1,length(tr));
    ar(1:2:end,:) = real(vr);  ar(2:2:end,:) = imag(vr);
    av = v'*ar;
    plot3(av(1:3:end-2,:)',av(2:3:end-1,:)',av(3:3:end,:)','.-'); hold all;
    grid on; pause; end

%    figure(1); clf; plot(tr,sum(abs(vr)-repmat(abs(veq),1,length(tr))),'.-'); pause;
%    figure(2); clf; plot(tr,ar(1:6,:),'.-'); grid on; pause; end

%% Plot RPOs in T vs Delta (include repeating orbits)
  clear;  load kse22orbits;  tmax = 200;  iskp = [112 119 155];
  trpo = [rpo.T];  phrpo = [rpo.d];  trpo(iskp) = [];  phrpo(iskp) = [];
  figure(1); clf; set(gcf,'pos',[420 280 465 660],'paperpos',[4.0 6.0 12.3 17.5]);
  plot([-phrpo phrpo], [trpo trpo], 'ro', [0 0], [0 tmax], 'k-');
  axis([-L/2 L/2 0 tmax]); hold on;
  for ii = 2:floor(tmax/trpo(1)),
    rt = trpo.*ii;  rph = mod(phrpo.*ii+L/2,L)-L/2; 
    irm = find(rt > tmax);  rt(irm) = [];  rph(irm) = [];
    plot([-rph rph], [rt rt], 'kx'); end
  print(gcf,'-dpng','ks22TvsDrep.png');
