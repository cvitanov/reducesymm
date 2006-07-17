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
  N = 63;  nu = (2*pi/54.0)^2;  h = 0.002;  randn('seed',12340000);
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
end

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
  
  d = 22;  N = 64;  x = (1:N)'./N*d;  h = 0.25;  tpre = 10;  t0 = 50;
%  a0 = zeros(N-2,1);  a0(5) = 1; %randn('seed',12340019);  a0(1:8) = 0.3*randn(8,1);
  [tt,aa] = ksfmedt(d, tpre, a0, h, 1);
  nap = zeros(size(aa,2),1);
  for ii = 1:size(aa,2), nap(ii) = norm(ksfm(t0,aa(:,ii),d)); end,
  [map,iap] = min(nap); tpre = tt(iap);
  figure(1); clf; plot(tt,nap,'.-',tpre,map,'ro'); pause;
  [tt,aa] = ksfmedt(d, tpre, a0, h);
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
    clear;  d = 12.5;  N = 32;  h = 0.25;
    x = d.*(1:N)'./N;  np = 1;  tend = 200;
    a0 = zeros(N-2,1); randn('seed',12340000);  a0(1:8) = 0.2*randn(8,1);
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
    options = optimset('Display','iter','Jacobian','off','TolFun',1e-9);
    a0 = fsolve(@(a)ksfmtr(a,d),[aa(:,ic); c(mc(ic))],options);
    for d = 12.5:.5:22;  disp(['L = ' num2str(d)]);
      a0 = fsolve(@(a)ksfmtr(a,d),a0,options);  end,
  
%%  Stability of the steady states and traveling wave
  clear; 
  load ks22uss3b; c = 0;
%   load ks22uss3a; c = 0;
%  c = -c;  a0(1:2:end) = -a0(1:2:end);
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
  clear;  d = 22;  N = 32;  h = 0.25;
  switch 2,  %%% Switch for seed selection
  case 1,  %%% Find close returns with a shift
    x = d.*(1:N)'./N;  np = 6;  tpre = 500;  tb = 10;  te = 150;
    if 0,
      a0 = zeros(N-2,1); randn('seed',12340002);  a0(1:8) = 0.2*randn(8,1);
    else
%      a0 = 0.1*ones(N-2,1); end,
%      load ks22utw1a; a0(1:2:end) = -a0(1:2:end);
      load ks22uqo074a;
%      load ks22seed47a;
    end
    [tt, aa, a0p, ap] = ksfmedt(d, tpre, a0, h, np);  na = size(aa,2);
    v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
    v = [zeros(1,size(aa,2)); v; zeros(1,size(aa,2)); flipud(conj(v))];
    u = real(fft(v)); 
    figure(1); set(gcf,'pos',[5 35 250 910]); clf; 
    pcolor(x,tt,u'); shading flat; tall(1.1); hold on;
    hp1 = plot(x([1 end]),tt([1 1]),'k-');
    figure(2); set(gcf,'pos',[835 525 560 420]); clf; 
    for ib = 0:1:200,
      nb = min(round(tb/(h*np))+ib, length(tt));  
      ne = min(round(te/(h*np))+nb+1, length(tt));
      set(hp1,'xdata',[x([1 end]);NaN;x([end 1])],'ydata',[tt([nb nb]) NaN tt([ne ne])]);
      a0 = aa(:,nb);  phase = (-0.5:0.01:0.49)'*d;
      v = aa(1:2:end,nb:ne)+1i*aa(2:2:end,nb:ne);  k = (1:N/2-1)';
      fac = exp((-2i*pi./d).*k*phase');  cr = zeros(size(fac,2),size(v,2));
      for ii = 1:size(fac,2),
        cr(ii,:) = sqrt(sum(abs(v.*repmat(fac(:,ii),1,size(v,2)) - ...
          repmat(v(:,1),1,size(v,2))).^2)); end
      contourf(tt(nb:ne)-tt(nb),phase,cr,[0:.1:1]); shading flat;
      [tend, ph, but] = ginput(1);  disp([ib but]); 
      if but == 32, break; end, end
  case 2,  %%% Load seed from a file
    load ks22seed115b;%  tend = 68;
%    load ks22uqo108a;
%    load ks22uss0a;
  end;
  switch 2, %%% Switch for method selection
  case 1,   %%% Find UQOs by associated flow
    alpha = [1e3 1e3];
    for ic = 1:20,
      switch 1,
      case 0,  C = eye(N-2);
      case 1, 
        [tt, aa, da] = ksfmjedt(d, tend, a0, h);
        ek = exp((-2i*pi/d).*ph.*(1:N/2-1)');
        ap = ek.*(aa(1:2:end,end)+1i*aa(2:2:end,end));
        disp(norm(ap-(a0(1:2:end)+1i*a0(2:2:end))));
        cda = repmat(ek,1,N-2).*(da(1:2:end,:)+1i*da(2:2:end,:));
        da(1:2:end,:) = real(cda); da(2:2:end,:) = imag(cda); 
        [u,C] = schur(da); disp(diag(C(1:6,1:6))');
        [u,C,v] = svd(da-eye(N-2));  C = -v*u';
      end,
      options = odeset('abstol',1e-7,'reltol',1e-8,'outputfcn',@ksdupoplot2);
      [s,y] = ode15s(@(t,a)ksfmf5new(t,a,d,h,C,alpha), [0 100000], ...
                     [a0; tend; ph], options);
      a0 = y(end-1,1:end-2)'; tend = y(end-1,end-1), ph = y(end-1,end),
      but = input('Stop? ','s');  if but == 'y', return; end, end
  case 2,   %%% Find UQOs using fsolve
    C = eye(N-2);  alpha = [1 1];  te = 0;
    options = optimset('Display','off','outputfcn',@ksoutfun,'MaxIter',2000, ...
              'MaxFunEvals',20000,'Jacobian','on','NonlEqnAlgorithm','lm','TolFun',1e-10);
    [upo,fval,eflag] = fsolve(@(a)ksfmf5new(te,a,d,h,C,alpha),[a0;tend;ph],options);
    a0 = upo(1:end-2); tend = upo(end-1), ph = upo(end),
  case 3,   %%% Find UQOs using lsqnonlin
    C = eye(N-2);  alpha = [1 1];
    options = optimset('Display','on','outputfcn',@ksoutfun,'MaxIter',2000,'Jacobian','on');
    lb = [repmat(-5,N-2,1);0;-d];  ub = [repmat(5,N-2,1);300;d];
    upo = lsqnonlin(@(a)ksfmf5new(te,a,d,h,C,alpha),[a0;tend;ph],lb,ub,options);
  end


%% List all located UPOs and UQOs
  clear;  lst = dir('ks22uqo10*.mat');  nst = length(lst);
  d = 22;  N = 32;  h = 0.25;  C = eye(N-2);  alpha = [1 1];
  for ist = 1:nst,
    load(lst(ist).name);  
    g = ksfmf5new(tend, [a0; tend; ph], d, h, C, alpha);
    disp(sprintf(' ''%15s'';%% %10.6f %10.6f %12.3e',lst(ist).name, tend, ph, norm(g(1:end-2))));
%    disp(sprintf(' ''%15s'';',lst(ist).name));
  end


%% Plot all located UPOs and UQOs
  clear;  name = [
  'ks22uqo016a.mat';%  16.316170   2.863787   1.726e-014
  'ks22uqo020a.mat';%  20.506740   0.000000   1.844e-007
  'ks22uqo033b.mat';%  32.806180  10.956763   9.982e-014
  'ks22uqo034a.mat';%  33.501034   4.045869   2.419e-013
  'ks22uqo035a.mat';%  34.639120   9.599483   1.825e-014
  'ks22uqo036a.mat';%  36.016825   3.911040   3.362e-011  <--
  
  'ks22uqo036d.mat';%  36.205721   3.813489   1.018e-013  <-- same orbit?
  'ks22uqo040b.mat';%  39.770736  -1.610607   1.812e-008
  'ks22uqo042a.mat';%  41.507477   6.479327   2.200e-008
  'ks22uqo046b.mat';%  46.052353   7.624041   2.406e-013  
  'ks22uqo047c.mat';%  46.502291  -7.758333   1.318e-013
  'ks22uqo047a.mat';%  47.323866   0.339464   8.331e-012
  
  'ks22uqo048a.mat';%  47.639322   5.675918   5.407e-010
  'ks22uqo056a.mat';%  55.599640  -5.247831   4.996e-014
  'ks22uqo060a.mat';%  59.892457  -5.442530   2.973e-012
  'ks22uqo066a.mat';%  65.612385   0.086473   1.671e-014  (33b x 2)
  'ks22uqo067a.mat';%  66.779323   0.000000   1.666e-006
  'ks22uqo068a.mat';%  67.939003   5.566339   6.073e-014
  
  'ks22uqo068e.mat';%  68.223923   1.706300   2.994e-008
  'ks22uqo069a.mat';%  68.537305   4.925039   4.496e-014
  'ks22uqo069b.mat';%  69.065162  -9.700812   3.349e-013
  'ks22uqo072a.mat';%  71.678083   5.502827   1.717e-014
  'ks22uqo074a.mat';%  73.980562   0.678530   8.188e-013
  'ks22uqo074c.mat';%  74.380749   4.434083   6.562e-009
  
  'ks22uqo075a.mat';%  75.483577  -6.235556   2.078e-009  
  'ks22uqo077a.mat';%  76.556365  -3.383947   3.993e-009  
  'ks22uqo077b.mat';%  76.604510  -1.677120   5.781e-010
  'ks22uqo078a.mat';%  78.194219   4.842480   3.825e-008
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
  'ks22uqo086b.mat';%  86.230815 -10.642810   1.816e-008
  'ks22uqo087a.mat';%  87.236323   0.000000   9.161e-007
  'ks22uqo089a.mat';%  89.419338  -0.788369   1.031e-007
  'ks22uqo091a.mat';%  91.354619   0.083919   1.993e-008
  
  'ks22uqo095a.mat';%  95.252982  -0.000001   1.432e-013
  'ks22uqo096b.mat';%  96.053538  -3.774012   5.968e-013
  'ks22uqo096a.mat';%  96.421035 -10.911410   9.371e-009  
  'ks22uqo098a.mat';%  98.254857  10.847250   4.403e-012
  'ks22uqo098b.mat';%  98.446990  -1.814935   6.338e-012
  'ks22uqo102a.mat';% 102.789193  -9.206446   1.026e-008

  'ks22uqo103a.mat';% 103.013668   1.006582   3.270e-013
  'ks22uqo108a.mat';% 108.078928 -10.270938   1.988e-011
  'ks22uqo109b.mat';% 109.115133  -6.089034   4.053e-008
  'ks22uqo109d.mat';% 109.288714   0.184272   3.416e-010
  'ks22uqo110a.mat';% 109.961306  -5.711897   1.547e-007  
  'ks22uqo110b.mat';% 110.253208   0.790874   1.053e-009
  
  'ks22uqo113a.mat';% 113.664192  -9.847374   8.846e-014
  'ks22uqo117a.mat';% 116.905523   6.640941   5.416e-009  
  'ks22uqo118a.mat';% 117.895472  -7.474551   5.115e-009
  'ks22uqo125a.mat';% 125.021554   9.062346   1.270e-009
  'ks22uqo127a.mat';% 127.126195   7.157910   9.227e-013
  'ks22uqo132b.mat';% 131.981245   1.342006   1.397e-012
  
  'ks22uqo133a.mat';% 132.696432  -0.015114   7.568e-013
  'ks22uqo139a.mat';% 139.219781   5.790933   4.777e-012  
  'ks22uqo140a.mat';% 140.031772   5.136071   5.650e-009
  'ks22uqo144a.mat';% 144.306458  -4.708727   2.288e-013  
  'ks22uqo171a.mat';% 170.936876   6.562272   1.063e-009
 ];  save ks22names name;
  nst = size(name,1);
  d = 22;  N = 32;  x = d.*(-N:N)'./(2*N);  h = 0.25;  np = 2;  te = 200;
  for ist = 1:nst,  ost = 1:nst;
%  for ist = 1:6, ost = [13:15 22 35 38];
    switch 1,
    case 1,
      ifig = floor(ist/6)+1;  iorb = mod(ist-1,6)+1;
      if iorb == 1,  figure(ifig); clf; orient landscape; 
        set(gcf,'pos',[35 200 1100 740]);
        hax = subplots(1,6,[0.03 0.03 0.07 0.1],[0.015 0 0 0]); end
      load(name(ost(ist),:));  ne = ceil(te./tend);
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
      axes(hax(iorb));  pcolor(x,tt,u'); shading flat; caxis([-3 3]);  hold on;
      plot(x([1 end])*ones(1,ne),[tend;tend]*(1:ne),'w-');
      plot(mod([ph;ph]*(1:ne-1),d)-d/2,[(1:ne-1);(2:ne)]*tend,'w-');
      title([sprintf('%3d: T=%7.3f ',ost(ist),tend) '\Delta=' sprintf('%9.6f',ph)]);
  %    title(sprintf('d/T = %8.4f',ph./tend));
      set(gca,'ticklength',0.7*get(gca,'ticklength'),'tickdir','out',...
        'ylim',[0 te],'ytick',20:20:te,'box','off');  xlabel('x');
      if iorb > 1, set(gca,'yticklabel',[]); end
      if iorb == 6,
        hax(7) = axes('pos',[.06 .9 .9 .07]);
        tit1 = 'Kuramoto-Sivashinsky: u_t = -uu_x - u_{xx} - u_{xxxx};  ';
        tit2 = 'x \in [0, L];    BC: u(x+L,t) = u(x,t);   L = 22.0;    ';
        tit3 = 'Solutions of the form: u(x+\Delta,T) = u(x,0)';
        text(0.1,.9,[tit1 tit2 tit3]);  set(hax(7),'visi','off');
        axes(hax(1)); ylabel('t','rotat',0); end
    case 2,  % Plot FMs of RPOs
      load(name(ist,:));
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

%% Introduce intermediate points (ksfmms0)
  clear; load ks22seed140a;  ni = [3 0];  d = 22;  N = 32;  h = 0.25;
  ni(2) = ceil(tend./(h.*(ni(1)+1)));
  [tt, aa] = ksfmedt(d, tend, a0, h, ni(2));
  ai = aa(:,1:end-1);  ai = [ai(:); tend-ni(1)*ni(2)*h; ph];
  gi = ksfmfns(ai, d, h, ni);
%  load ks22seed123ai;
  options = optimset('Display','on', 'outputfcn',@ksoutfun,'MaxIter',2000, ...
      'MaxFunEvals',100000,'Jacobian','off','NonlEqnAlgorithm','lm','TolFun',1e-9);
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