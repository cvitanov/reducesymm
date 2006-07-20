%% Single movie plot (1st attempt at making avi movie)
if 0,
  load ks092mov;  ip = 1:1:size(y,1);
  hfig = figure(1); set(gcf,'pos',[30   450   530   500]); clf;
%  set(hfig,'DoubleBuffer','on');
%  mov = avifile('ks092mov.avi','fps',10,'quality',100);
  [tt, ye] = kscvfm(nu, y(1,end),y(1,1:end-1)',h,2);
  ng = norm(ye(:,end)-ye(:,1));
  hp2 = plot(ye(1,:),ye(2,:),'-b.'); hold on;
  hp1 = plot(y(1,1),y(1,2),'ro','markersize',4);
  hp3 = plot(ye(1,end),ye(2,end),'ko','markersize',4);
  hp4 = text(0.6,-1.3,sprintf('|g| = %8.6f',ng));
  hp5 = text(0.6,-1.5,sprintf(' T = %8.6f',y(1,end)));
  axis image; axis([-1.2 1.3 -1.7 0.8]);
  set(gca,'pos',[0.11 0.12 0.78 0.81],'xtick',-1:1,'ytick',-1:0);
  xlabel('a_1','fontsize',14,'fontname','times','fontangle','italic'); 
  ylabel('a_2','fontsize',14,'fontname','times','fontangle','italic','rotat',0);
  title('Kuramoto-Sivashinsky (FM):  \nu = 0.015,  d = 31','fontsize',12,'fontname','times');
%  hax = axes('pos',[0 0 1 1]); set(hax,'visible','off'); pause;
  j = 0;
  for i = ip,
    [tt, ye] = kscvfm(nu, y(i,end),y(i,1:end-1)',h,2);
    ng = norm(ye(:,end)-ye(:,1));
    set(hp2,'xdata',ye(1,:),'ydata',ye(2,:));
    set(hp1,'xdata',[get(hp1,'xdata') y(i,1)],'ydata',[get(hp1,'ydata') y(i,2)]);
    set(hp3,'xdata',[get(hp3,'xdata') ye(1,end)],'ydata',[get(hp3,'ydata') ye(2,end)]);
    set(hp4,'string',sprintf('|g| = %8.6f',ng));
    set(hp5,'string',sprintf(' T = %8.6f',y(i,end)));
    j = j+1;  M(j) = getframe(gcf);
%    F = getframe(hax);  mov = addframe(mov,F);
  end,
%  mov = close(mov);
end,

%% Movie for UPOs (odd solutions, ks092 and ks099)
if 0, clear;  name = 'ks099';
  load([name 'mov']);  log = load([name '.dat']);
  log(:,1) = log(:,1)-1;  log(:,3) = log(:,3).*5;  log(:,4) = log(:,4)./nu; %%% scale
  figure(1); clf;  set(gcf,'pos',[100 300 600 632]);
  hax1 = axes('pos',[0.08 0.14 0.36 0.82]);
  hax2 = axes('pos',[0.53 0.14 0.43 0.36]);
  hax3 = axes('pos',[0.53 0.60 0.43 0.36]);
  [tt, ye] = kscvfm(nu, y(1,end),y(1,1:end-1)',h,3);  tt = tt./nu; %%% scale
  Nu = 128;  N = size(ye,1);  na = size(ye,2);
  ue = real(ifft(1i*(2*N+2)*[zeros(1,na); ye; zeros(Nu-1-2*N,na); -flipud(ye)]));
  ue = ue(1:Nu/2+1,:);  xx = 2*pi*(0:Nu/2)'/Nu;  xx = xx./sqrt(nu); %%% scale
 axes(hax1); pcolor(xx,tt,ue'); shading interp;  
  xlabel('x','fontsize',14,'fontname','times'); 
  ylabel('t','fontsize',14,'fontname','times','rotat',0);
 axes(hax2); 
  ng = sqrt(sum((ye-repmat(ye(:,1),1,size(ye,2))).^2));
  hp21 = plot(tt,ng,'-g.','color',[0 .7 0]); hold on;
  hp22 = plot(tt(end),ng(end),'ko','markersize',5);
  hp23 = text(-80,-1.8,sprintf('step = %3d   nfe = %3d   s = %7.3f   T = %7.4f  |g| = %8.6f',log(1,:)),'fontname','courier');
  axis([0 73 0 6.5]);  set(gca,'ytick',0:2:6);
  xlabel('t','fontsize',14,'fontname','times'); ylabel('|a(t) - a(0)|','fontsize',14,'fontname','times');
 axes(hax3);
  hp31 = plot(ye(1,:),ye(2,:),'-b.'); hold on;
  hp32 = plot(y(1,1),y(1,2),'ro','markersize',5);
  hp33 = plot(ye(1,end),ye(2,end),'ko','markersize',5);
  axis([-1.2 1.3 -1.7 0.8]);  set(gca,'xtick',-1:1,'ytick',-1:0);
  xlabel('a_1','fontsize',14,'fontname','times');  ylabel('a_2','fontsize',14,'fontname','times','rotat',0);
  jj = 1;  clear M;  M(jj) = getframe(gcf);
%  for ii = 2:size(y,1), % ks092
  for ii = [2:80 82:3:222 227:5:size(y,1)],  % ks099
    [tt, ye] = kscvfm(nu, y(ii,end),y(ii,1:end-1)',h,3);  tt = tt./nu; %%% scale
    Nu = 128;  N = size(ye,1);  na = size(ye,2);
    ue = real(ifft(1i*(2*N+2)*[zeros(1,na); ye; zeros(Nu-1-2*N,na); -flipud(ye)]));
    ue = ue(1:Nu/2+1,:);%  xx = 2*pi*(0:Nu/2)'/Nu;
   axes(hax1); pcolor(xx,tt,ue'); shading interp;
    xlabel('x','fontsize',14,'fontname','times'); 
    ylabel('t','fontsize',14,'fontname','times','rotat',0);
   axes(hax2);
    ng = sqrt(sum((ye-repmat(ye(:,1),1,size(ye,2))).^2));
    set(hp21,'xdata',tt,'ydata',ng);  
    set(hp22,'xdata',[get(hp22,'xdata') tt(end)],'ydata',[get(hp22,'ydata') ng(end)]);
    set(hp23,'string',sprintf('step = %3d   nfe = %3d   s = %7.3f   T = %7.4f  |g| = %8.6f',log(ii,:)));
   axes(hax3);
    set(hp31,'xdata',ye(1,:),'ydata',ye(2,:));
    set(hp32,'xdata',[get(hp32,'xdata') ye(1,1)],'ydata',[get(hp32,'ydata') ye(2,1)]);
    set(hp33,'xdata',[get(hp33,'xdata') ye(1,end)],'ydata',[get(hp33,'ydata') ye(2,end)]);
    jj = jj + 1;  M(jj) = getframe(gcf);
  end,
%  movie2avi(M,[name 'mov.avi'],'fps',5,'quality',100);
end,

%% Movie for RPO detection
if 1,  clear;  name = 'ks22uqo042';  N = 32;  np = 2;
  load([name 'mov']);  log = load([name '.dat']);  x = d.*(-N:N)'./(2*N);
  figure(1); clf;  set(gcf,'pos',[200 300 600 632]);
  hax1 = axes('pos',[0.08 0.14 0.35 0.82]);
  hax2 = axes('pos',[0.53 0.14 0.43 0.36]);
  hax3 = axes('pos',[0.53 0.60 0.43 0.36]);
  [tt, aa] = ksfmedt(d, y(1,N-1), y(1,1:N-2)', h, np);  nt = length(tt);
  v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
  vv = [zeros(1,nt); v; zeros(N+1,nt); flipud(conj(v))];
  u = real(fft(vv));  u = [u; u(1,:)];
  ek = exp((-2i*pi*y(1,N)./d).*(1:N/2-1)');  vp = v(:,end).*ek;
  ap = zeros(N-2,1);  ap(1:2:end) = real(vp);  ap(2:2:end) = imag(vp);
  ng = sqrt(sum(([aa ap]-repmat(aa(:,1),1,nt+1)).^2));
 axes(hax1); pcolor(x,tt,u'); shading interp;  caxis([-3 3]);
  xlabel('x','fontsize',14,'fontname','times');
  ylabel('t','fontsize',14,'fontname','times','rotat',0);
 axes(hax2); 
  hp21 = plot([tt tt(end)],ng,'-g.','color',[0 .7 0]); hold on;
  hp22 = plot(tt(end),ng(end),'ko','markersize',5);
  hp23 = text(-50,-0.4,[sprintf('step = %3d;  s = %7.3f;  T = %7.4f;',log(1,[1 3 4]))...
        '  \Delta = ' sprintf('%6.4f;  |g| = %8.6f',log(1,5:6))],'fontname','courier');
  axis([0 45 0 1.5]);  set(gca,'ytick',0:.5:1);
  xlabel('t','fontsize',14,'fontname','times'); 
  ylabel('|u(t) - u(0)|','fontsize',14,'fontname','times');
 axes(hax3);
  hp31 = plot([aa(1,:) ap(1)],[aa(2,:) ap(2)],'-b.'); hold on;
  hp32 = plot(aa(1,1),aa(2,1),'ro','markersize',5);
  hp33 = plot(ap(1),ap(2),'ko','markersize',5);
  axis([-0.55 0.3 -0.3 0.2]);  set(gca,'xtick',-.4:.2:.2,'ytick',-.2:.1:.1);
  xlabel('Re u_1','fontsize',14,'fontname','times');  
  ylabel('Im u_1','fontsize',14,'fontname','times','verticalalign','baseline');
  jj = 1;  clear M;  M(jj) = getframe(gcf);
  for ii = [2:2:80 81:120 122:2:size(y,1)],
    [tt, aa] = ksfmedt(d, y(ii,N-1), y(ii,1:N-2)', h, np);  nt = length(tt);
    v = aa(1:2:end,:) + 1i*aa(2:2:end,:);
    vv = [zeros(1,nt); v; zeros(N+1,nt); flipud(conj(v))];
    u = real(fft(vv));  u = [u; u(1,:)];
    ek = exp((-2i*pi*y(ii,N)./d).*(1:N/2-1)');  vp = v(:,end).*ek;
    ap = zeros(N-2,1);  ap(1:2:end) = real(vp);  ap(2:2:end) = imag(vp);
    ng = sqrt(sum(([aa ap]-repmat(aa(:,1),1,nt+1)).^2));
   axes(hax1); pcolor(x,tt,u'); shading interp; caxis([-3 3]);
    xlabel('x','fontsize',14,'fontname','times'); 
    ylabel('t','fontsize',14,'fontname','times','rotat',0);
   axes(hax2);
    set(hp21,'xdata',[tt tt(end)],'ydata',ng);  
    set(hp22,'xdata',[get(hp22,'xdata') tt(end)],'ydata',[get(hp22,'ydata') ng(end)]);
    set(hp23,'string',[sprintf('step = %3d;  s = %7.3f;  T = %7.4f;',log(ii,[1 3 4]))...
        '  \Delta = ' sprintf('%6.4f;  |g| = %8.6f',log(ii,5:6))]);
   axes(hax3);
    set(hp31,'xdata',[aa(1,:) ap(1)],'ydata',[aa(2,:) ap(2)]);
    set(hp32,'xdata',[get(hp32,'xdata') aa(1,1)],'ydata',[get(hp32,'ydata') aa(2,1)]);
    set(hp33,'xdata',[get(hp33,'xdata') ap(1)],'ydata',[get(hp33,'ydata') ap(2)]);
    jj = jj + 1;  M(jj) = getframe(gcf);
  end,
  movie2avi(M,[name 'mov.avi'],'fps',5,'quality',100);
end,