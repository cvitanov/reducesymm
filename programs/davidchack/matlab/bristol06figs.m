%% Long time series (L = 50)
  clear;  d = 50;  N = 64;  h = 0.25;  x = d.*(-N:N)'./(2*N);
  a0 = zeros(N-2,1); randn('seed',12340001);  a0(1:12) = 0.2*randn(12,1);
  tpre = 50;  tend = 600;  np = 2;  nax = 3;
  [tt, aa] = ksfmedt(d, tpre, a0, h);  a0 = aa;  % pre-iterates
  [tt, aa] = ksfmedt(d, tend, a0, h, np);  nt = length(tt);
  v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
  vv = [zeros(1,nt); v; zeros(N+1,nt); flipud(conj(v))];
  u = real(fft(vv));  u = [u; u(1,:)];
  ik = -(2i*pi/d).*(1:N/2-1)';  vx = repmat(ik,1,nt).*v;
  vv = [zeros(1,nt); vx; zeros(N+1,nt); flipud(conj(vx))];
  ux = real(fft(vv)); ux = [ux; ux(1,:)];
  fig1 = figure('PaperOrientation','landscape',...
    'PaperPosition',[0.6345 0.6345 28.41 19.72],...
    'PaperSize',[29.68 20.98],'Position',[200  270  900  400]);
  hax = subplots(1,nax,[0.01 0.03 0.12 0.03],[0.05 0 0 0]);
  for ia = 1:nax,  
    it = ia*floor(tend/h)/(nax*np) + (-floor(tend/h)/(nax*np):0)+1;
    axes(hax(ia));
    pcolor(x,tt(it),u(:,it)'); shading interp;  caxis([-3 3]);
    xlabel('$x$','fontsize',15,'interp','latex'); end
  axes(hax(1));  ylabel('Time','fontsize',14);
  
%% Long time series (antisymmetric subspace)
  clear;  d = 51.3;  N = 64;  h = 0.25;  x = d.*(-N:N)'./(2*N);
  a0 = zeros(N-2,1); randn('seed',12340002);  a0(2:2:8) = 0.2*randn(4,1);
  tpre = 50;  tend = 1000;  np = 2;  nax = 5;
  [tt, aa] = ksfmedt(d, tpre, a0, h);  a0 = aa;  % pre-iterates
  [tt, aa] = ksfmedt(d, tend, a0, h, np);  nt = length(tt);
  v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
  vv = [zeros(1,nt); v; zeros(N+1,nt); flipud(conj(v))];
  u = real(fft(vv));  u = [u; u(1,:)];
  ik = -(2i*pi/d).*(1:N/2-1)';  vx = repmat(ik,1,nt).*v;
  vv = [zeros(1,nt); vx; zeros(N+1,nt); flipud(conj(vx))];
  ux = real(fft(vv)); ux = [ux; ux(1,:)];
  fig1 = figure('PaperOrientation','landscape',...
    'PaperPosition',[0.6345 0.6345 28.41 19.72],...
    'PaperSize',[29.68 20.98],'Position',[200  270  900  400]);
  hax = subplots(1,nax,[0.01 0.03 0.12 0.03],[0.05 0 0 0]);
  for ia = 1:nax,  
    it = ia*floor(tend/h)/(nax*np) + (-floor(tend/h)/(nax*np):0)+1;
    axes(hax(ia));
    pcolor(x,tt(it),u(:,it)'); shading interp;  caxis([-3 3]);
    xlabel('$x$','fontsize',15,'interp','latex'); end
  axes(hax(1));  ylabel('Time','fontsize',14);
  
%% Long time series (full space)
  clear;  d = 22;  N = 32;  h = 0.25;  x = d.*(-N:N)'./(2*N);
  a0 = zeros(N-2,1); randn('seed',12340000);  a0(1:8) = 0.2*randn(8,1);
  tpre = 50;  tend = 1000;  np = 2;  nax = 5;
  [tt, aa] = ksfmedt(d, tpre, a0, h);  a0 = aa;  % pre-iterates
  [tt, aa] = ksfmedt(d, tend, a0, h, np);  nt = length(tt);
  v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
  vv = [zeros(1,nt); v; zeros(N+1,nt); flipud(conj(v))];
  u = real(fft(vv));  u = [u; u(1,:)];
  ik = -(2i*pi/d).*(1:N/2-1)';  vx = repmat(ik,1,nt).*v;
  vv = [zeros(1,nt); vx; zeros(N+1,nt); flipud(conj(vx))];
  ux = real(fft(vv)); ux = [ux; ux(1,:)];
  fig1 = figure('PaperOrientation','landscape',...
  'PaperPosition',[0.6345 0.6345 28.41 19.72],...
  'PaperSize',[29.68 20.98],'Position',[200  270  1000  680]);
  for ia = 1:nax,  
    it = ia*floor(tend/h)/(nax*np) + (-floor(tend/h)/(nax*np):0)+1;
    ax(ia) = axes('position',[0.19*ia-0.13 0.09 0.15 0.85]);
    pcolor(x,tt(it),u(:,it)'); shading interp;  caxis([-3 3]);
    xlabel('x','fontsize',14); end
  axes(ax(1));  ylabel('Time','fontsize',14);
  
%% Steady states and traveling waves
  clear;  d = 22;  N = 32;  h = 0.25;  x = d.*(-N:N)'./(2*N);
  fig1 = figure('PaperOrientation','landscape',...
  'PaperPosition',[0.6345 0.6345 28.41 19.72],...
  'PaperSize',[29.68 20.98],'Position',[200  270  800  680]);
  for ia = 1:4,
    ax(ia) = axes('position',[0.24*ia-0.17 0.09 0.19 0.85]);  end
 axes(ax(1));  load ks22uss0a;   %  Steady 2-wave state
   if 1, tend = 200;  np = 4;
     [tt, aa] = ksfmedt(d, tend, a0, h, np);  nt = length(tt);
     v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
     vv = [zeros(1,nt); v; zeros(N+1,nt); flipud(conj(v))];
     u = real(fft(vv));  u = [u; u(1,:)];
   else
     aa = [a0 a0];  tt = [0 200]; end
   pcolor(x,tt,u'); shading interp; caxis([-3 3]);
   xlabel('x','fontsize',14);  ylabel('Time','fontsize',14);
 axes(ax(2));  load ks22uss3a;   %  Steady 3-wave state
   if 1, tend = 200;  np = 4;
     [tt, aa] = ksfmedt(d, tend, a0, h, np);  nt = length(tt);
     v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
     vv = [zeros(1,nt); v; zeros(N+1,nt); flipud(conj(v))];
     u = real(fft(vv));  u = [u; u(1,:)];
   else
     aa = [a0 a0];  tt = [0 200]; end
   pcolor(x,tt,u'); shading interp; caxis([-3 3]);
   xlabel('x','fontsize',14);
 axes(ax(3));  load ks22utw1a;   %  Traveling 1-wave
   if 1, tend = 200;  np = 1;
     [tt, aa] = ksfmedt(d, tend, a0, h, np);  nt = length(tt);
     v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
     vv = [zeros(1,nt); v; zeros(N+1,nt); flipud(conj(v))];
     u = real(fft(vv));  u = [u; u(1,:)];
   else
     tt = 0:0.5:200;  nt = length(tt);
     v = repmat(a0(1:2:end)+1i*a0(2:2:end),1,nt);
     v = exp(2i*pi/d*c*(1:15)'*tt).*v;
     vv = [zeros(1,nt); v; zeros(N+1,nt); flipud(conj(v))];
     u = real(fft(vv));  u = [u; u(1,:)];   end
   pcolor(x,tt,u'); shading interp; caxis([-3 3]);
   xlabel('x','fontsize',14);
 axes(ax(4));  load ks22utw1a;   %  Traveling 1-wave (reverse direction)
   if 1, tend = 200;  np = 1;  a0(1:2:end) = -a0(1:2:end);
     [tt, aa] = ksfmedt(d, tend, a0, h, np);  nt = length(tt);
     v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
     vv = [zeros(1,nt); v; zeros(N+1,nt); flipud(conj(v))];
     u = real(fft(vv));  u = [u; u(1,:)];
   else
     a0(1:2:end) = -a0(1:2:end);  c = -c;
     tt = 0:0.5:200;  nt = length(tt);
     v = repmat(a0(1:2:end)+1i*a0(2:2:end),1,nt);
     v = exp(2i*pi/d*c*(1:15)'*tt).*v;
     vv = [zeros(1,nt); v; zeros(N+1,nt); flipud(conj(v))];
     u = real(fft(vv));  u = [u; u(1,:)];   end
   pcolor(x,tt,u'); shading interp;  caxis([-3 3]);
   xlabel('x','fontsize',14);

%%  Relative periodic orbits (plotted in ksdupo)

%% Plot FMs of RPOs in co-rotating frame
  clear;  d = 22;  N = 32;  x = d.*(-N:N)'./(2*N);  h = 0.25;  ne = 2;
  load ks22uqo060a;
  [tti, aai] = ksfmedt(d, tend, a0, h, 1);
  ek = exp((2i*pi/d).*ph.*(1:N/2-1)');
  tt = tti(1:end-1);  aa = aai(:,1:end-1);
  for ie = 1:ne-1,
    vi = (aai(1:2:end,:)+1i*aai(2:2:end,:)).*repmat(ek,1,size(aai,2));
    aai(1:2:end,:) = real(vi);  aai(2:2:end,:) = imag(vi);
    aa = [aa aai(:,1:end-1)];  tt = [tt tti(1:end-1)+ie*tend]; end
  v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
  vv = [zeros(1,size(aa,2)); v; zeros(N+1,size(aa,2)); flipud(conj(v))];
  u = real(fft(vv)); u = [u; u(1,:)];
fig1 = figure('PaperPosition',[0.6345 6.345 20.3 15.23],...
    'PaperSize',[20.98 29.68],'Position',[5  500  1046  440]);
 ax1 = axes('Position',[0.0526 0.0923 0.149 0.85],'Parent',fig1);
  pcolor(x,tt,u'); shading flat;  caxis([-3 3]);  hold on; 
  plot(x([1 end])*ones(1,ne),[tend;tend]*(1:ne),'w-');
  plot(mod([ph;ph]*(1:ne-1),d)-d/2,[(1:ne-1);(2:ne)]*tend,'w-');
  title([sprintf('T = %6.2f  ',tend) '\Delta = ' sprintf('%9.6f',ph)]);  
  xlabel('x'); ylabel('Time');

  ek = exp((-2i*pi/d).*ph./tend*(1:N/2-1)'*tt);  vr = ek.*v;
  e1 = exp(-2i*pi./tend*(1:N/2-1)'*tt);
  load ks22uss0a;  v2 = repmat(a0(1:2:end)+1i*a0(2:2:end),1,length(tt));
  load ks22uss3a;  v3 = repmat(a0(1:2:end)+1i*a0(2:2:end),1,length(tt));
  load ks22utw1a;  v4 = repmat(a0(1:2:end)+1i*a0(2:2:end),1,length(tt));
  v2r = e1.*v2;  v3r = e1.*v3;  v4r = e1.*v4;  v5r = -e1.*conj(v4);
 ax2 = axes('Position',[0.261 0.142 0.242 0.85],'Parent',fig1);
  axis(ax2,[-0.8119 0.5598 -0.5598 0.5598 -1.24 1.24]);
  plot3(real(vr(2,:)), imag(vr(2,:)), imag(vr(3,:)), '.-');
  grid on; hold on; 
  xlabel('Re u_2');  ylabel('Im u_2');  zlabel('Re u_3');
  plot3(real(v2r(2,:)), imag(v2r(2,:)), imag(v2r(3,:)), 'k.-');
  plot3(real(v3r(2,:)), imag(v3r(2,:)), real(v3r(3,:)), 'r.-'); 
  axis image;  view(-160,30);
 ax3 = axes('Position',[0.589 0.129 0.397 0.85],'Parent',fig1);
  axis(ax3,[-1.24 1.24 -1.24 1.24 -0.8119 0.5598]);
  plot3(real(vr(3,:)), imag(vr(3,:)), real(vr(2,:)), '.-');
  grid on; hold on;
  xlabel('Re u_3');  ylabel('Im u_3');  zlabel('Re u_2');
  plot3(real(v2r(3,:)), imag(v2r(3,:)), real(v2r(2,:)), 'k.-');
  plot3(real(v3r(3,:)), imag(v3r(3,:)), real(v3r(2,:)), 'r.-');
  axis image;  view(-160,30);

%% Plot FMs of RPOs in co-rotating frame
  clear;  d = 22;  N = 32;  x = d.*(-N:N)'./(2*N);  h = 0.25;  ne = 2;
  load ks22uqo077b;
  [tti, aai] = ksfmedt(d, tend, a0, h, 1);
  ek = exp((2i*pi/d).*ph.*(1:N/2-1)');
  tt = tti(1:end-1);  aa = aai(:,1:end-1);
  for ie = 1:ne-1,
    vi = (aai(1:2:end,:)+1i*aai(2:2:end,:)).*repmat(ek,1,size(aai,2));
    aai(1:2:end,:) = real(vi);  aai(2:2:end,:) = imag(vi);
    aa = [aa aai(:,1:end-1)];  tt = [tt tti(1:end-1)+ie*tend]; end
  v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
  vv = [zeros(1,size(aa,2)); v; zeros(N+1,size(aa,2)); flipud(conj(v))];
  u = real(fft(vv)); u = [u; u(1,:)];
fig1 = figure('PaperPosition',[0.6345 6.345 20.3 15.23],...
    'PaperSize',[20.98 29.68],'Position',[5  390  1046  550]);
 ax1 = axes('Position',[0.0526 0.0923 0.149 0.85],'Parent',fig1);
  pcolor(x,tt,u'); shading flat;  caxis([-3 3]);  hold on; 
  plot(x([1 end])*ones(1,ne),[tend;tend]*(1:ne),'w-');
  plot(mod([ph;ph]*(1:ne-1),d)-d/2,[(1:ne-1);(2:ne)]*tend,'w-');
  title([sprintf('T = %6.2f  ',tend) '\Delta = ' sprintf('%9.6f',ph)]);  
  set(gca,'ytick',0:20:ne*tend); xlabel('x'); ylabel('Time');

  ek = exp((-2i*pi/d).*ph./tend*(1:N/2-1)'*tt);  vr = ek.*v;
  e1 = exp(-2i*pi./tend*(1:N/2-1)'*tt);
  load ks22uss0a;  v2 = repmat(a0(1:2:end)+1i*a0(2:2:end),1,length(tt));
  load ks22uss3a;  v3 = repmat(a0(1:2:end)+1i*a0(2:2:end),1,length(tt));
  load ks22utw1a;  v4 = repmat(a0(1:2:end)+1i*a0(2:2:end),1,length(tt));
  v2r = e1.*v2;  v3r = e1.*v3;  v4r = e1.*v4;  v5r = -e1.*conj(v4);
 ax2 = axes('Position',[0.274 0.100 0.242 0.85],'Parent',fig1);
  axis(ax2,[-0.56 0.56 -0.826 0.56 -1.24 1.24]);
  plot3(real(vr(2,:)), imag(vr(2,:)), imag(vr(3,:)), '.-');
  grid on; hold on; 
  xlabel('Re u_2');  ylabel('Im u_2');  zlabel('Re u_3');
  plot3(real(v2r(2,:)), imag(v2r(2,:)), imag(v2r(3,:)), 'k.-');
  plot3(real(v3r(2,:)), imag(v3r(2,:)), real(v3r(3,:)), 'r.-'); 
  axis image;  view(-60,30);
 ax3 = axes('Position',[0.598 0.129 0.389 0.85],'Parent',fig1);
  axis(ax3,[-1.24 1.24 -1.24 1.24 -0.56 0.56]);
  plot3(real(vr(3,:)), imag(vr(3,:)), real(vr(2,:)), '.-');
  grid on; hold on;
  xlabel('Re u_3');  ylabel('Im u_3');  zlabel('Re u_2');
  plot3(real(v2r(3,:)), imag(v2r(3,:)), real(v2r(2,:)), 'k.-');
  plot3(real(v3r(3,:)), imag(v3r(3,:)), real(v3r(2,:)), 'r.-');
  axis image;  view(60,30);

%% Plot all detected RPOs as dots in T vs. phase  
  load kse22orbits;
  fig1 = figure('PaperPosition',[0.6345 6.345 20.3 15.23],...
         'PaperSize',[20.98 29.68],'Position',[5  300  350 550]);
  ax1 = axes('Position',[0.14 0.08 0.80 0.90],'Parent',fig1);
  plot([1;-1]*[rpo.d],[1;1]*[rpo.T],'ro','markersize',5); hold on;
  plot([0; 0],[0 200],'k-'); axis([-L/2 L/2 0 200]);
  