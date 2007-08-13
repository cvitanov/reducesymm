%% E1 1st unstable manifold and orbits on it
  clear;  load kse22orbits;  k = 1;  p = 1;  h = 0.1;  tend = 75;
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,6))]);
  for delta = [0:0.06:ere*period 0.3 1.2],
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmedt(L, tend-delta/ere, a0, h, 2);
    if delta == 0, av = v'*aa;
    else aa = [aa repmat(NaN,size(aa,1),size(av,2)-size(aa,2))]; av = [av; v'*aa]; end
    if delta == 0.3, tt1 = tt; aa1 = aa(:,1:length(tt1)); end, end
  figure(1); set(gcf,'pos',[5 400 420 300],'paperpos',[5 8 12 9]); clf;
  ax1 = axes('pos',[0.04 0.08 0.92 0.90]);
    plot3(av(1:3:end-8,1:2:end)',av(2:3:end-7,1:2:end)',av(3:3:end-6,1:2:end)','-','color',[.5 .5 .5]);
    hold on; grid on; axis equal; view(-30,35);
    plot3(av(end-5,:)',av(end-4,:)',av(end-3,:)','r.-');
    plot3(av(end-2,:)',av(end-1,:)',av(end,:)','.-','color',[0 .8 0]);
    xlabel('$v_1$','fontsize',14,'interp','latex','pos',[0 -0.9 -0.4]);
    ylabel('$v_2$','fontsize',14,'interp','latex','pos',[-0.9 0.0 -0.5]);
    zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex','pos',[-0.9 0.6 0]);
    text(0.7,0.0,0.2,'B','fontsize',14,'color','r');
    text(-0.5,0.5,0.1,'A','fontsize',14,'color',[0 .8 0]);
    text(-0.1, 0.0,-0.25,'$E_1$','fontsize',14,'interp','latex');    
    print -depsc2 ks22_E1_plane1_manifold_c.eps
    figure(2); clf; set(gcf,'pos',[450 400 420 300],'paperpos',[5 8 9 9]);
  ax2 = axes('pos',[0.58 0.13 0.34 0.75]);
    a0 = eq(k).a + 1e-4.*exp(0.3).*v(:,2);
    [tt, aa] = ksfmedt(L, tend, a0, h, 5);  [x, uu] = ksfm2real(aa(1:30,:), L); 
    pcolor(x,tt(1:2:end),uu(:,1:2:end)'); caxis([-3 3]); shading flat;
    xlabel('$x$','fontsize',14,'interp','latex');
    ht = title('Orbit B');  set(ht,'color','r','fontsize',15);
  ax3 = axes('pos',[0.13 0.13 0.34 0.75]);
    a0 = eq(k).a + 1e-4.*exp(1.2).*v(:,2);
    [tt, aa] = ksfmedt(L, tend, a0, h, 5);  [x, uu] = ksfm2real(aa(1:30,:), L);
    pcolor(x,tt(1:2:end),uu(:,1:2:end)'); caxis([-3 3]); shading flat;  
    ylabel('$t$','rotat',0,'fontsize',14,'interp','latex'); 
    xlabel('$x$','fontsize',14,'interp','latex');
    ht = title('Orbit A');  set(ht,'color',[0 .8 0],'fontsize',15);
    print -depsc2 ks22_E1_plane1_orbits_c.eps

%% E1 evolving into full chaos
  clear;  load kse22orbits;  tend = 150;  np = 10;
  a1 = eq(1).a + 1e-3.*real(eq(1).evec(:,1));
  [tt, aa] = ksfmedt(L, tend, a1, h, np);  [x, uu] = ksfm2real(aa, L, 64);
  figure(1); set(gcf,'pos',[100 500 250 400]); clf;
  ax1 = axes('pos',[0.22 0.15 0.70 0.80]); pcolor(x,tt,uu'); caxis([-3 3]);
  shading interp; colormap('jet');
  xlabel('x','fontsize',15);  ylabel('t','fontsize',15,'rotat',0);
  set(get(gca,'ylabel'),'pos',[-13.5 74 1]);
  set(gcf,'paperpos',[8 10 6 10]); 
%  print -depsc2 cvitanovic/kse22_E1_chaos_color.eps
  
%% E1 2nd unstable manifold
  clear;  load kse22orbits;  k = 1;  h = 0.1;  tend = 220;  av = [];
  ere = real(eq(k).eig(3));  period = 2*pi/imag(eq(k).eig(3));
  v = gsorth([real(eq(k).evec(:,3)) imag(eq(k).evec(:,3)) real(eq(k).evec(:,7))]);
  for delta = [0:0.05:ere*period 0.8 0.098 0.1028],
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmedt(L, tend, a0, h, 4); av = [av; v'*aa];
    if delta == 0.8, aa1 = aa; end
    if delta == 0.098, aa2 = aa; end, end,
  figure(1); set(gcf,'pos',[100 200 700 700]); clf;
  ax1 = axes('pos',[0.1 0.46 0.83 0.52]);
    plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','-','color',[.5 .5 .5]);
    hold on; grid on; axis equal; view(-110,30);
    plot3(av(end-8,:)',av(end-7,:)',av(end-6,:)','k-','linewidth',1.8); % blue
    plot3(av(end-5,:)',av(end-4,:)',av(end-3,:)','k:','linewidth',2.9); % red
    plot3(av(end-2,:)',av(end-1,:)',av(end,:)','k--','linewidth',1.8);  % green
    e2 = v'*aa(:,end); plot3(e2(1),e2(2),e2(3),'k.','markersize',15);
    e3 = v'*eq(3).a; plot3(e3(1),e3(2),e3(3),'k.','markersize',15);
    text(-0.4,-0.1,-0.2,'A','fontsize',15);
    text(-0.4,-0.5,-0.3,'B','fontsize',15);
    text(-0.4,-0.5, 0.6,'C','fontsize',15);
    text(-0.4, 0.3, 0.52,'E2','fontsize',15);    
    text(-0.4,-0.75, 0.35,'E3','fontsize',15);    
  ax2 = axes('pos',[0.10 0.08 0.23 0.28]);
    [x, uu] = ksfm2real(aa1(:,1:2:end), L); pcolor(x,tt(1:2:end),uu'); caxis([-3 3]);
    shading interp;  xlabel('x','fontsize',15);  ylabel('t','rotat',0,'fontsize',15);
    ht = title('Orbit A','fontsize',15);
  ax3 = axes('pos',[0.41 0.08 0.23 0.28]);
    [x, uu] = ksfm2real(aa2(:,1:2:end), L); pcolor(x,tt(1:2:end),uu'); caxis([-3 3]);
    shading interp;  xlabel('x','fontsize',15);%  ylabel('t','rotat',0,'fontsize',15);
    ht = title('Orbit B','fontsize',15);
  ax4 = axes('pos',[0.72 0.08 0.23 0.28]);
    [x, uu] = ksfm2real(aa(:,1:2:end), L);  pcolor(x,tt(1:2:end),uu'); caxis([-3 3]);
    shading interp;  xlabel('x','fontsize',15);%  ylabel('t','rotat',0,'fontsize',15);
    ht = title('Orbit C','fontsize',15);
    set(gcf,'paperpos',[3 6 14 17]); wysiwyg;
%    print -depsc2 cvitanovic/kse22_E1_UM2.eps
%    print -dpdf  cvitanovic/kse22_E1_UM2.pdf

%% E1 2nd unstable manifold (without orbits)
  clear;  load kse22orbits;  k = 1;  h = 0.1;  tend = 220;  av = [];
  ere = real(eq(k).eig(3));  period = 2*pi/imag(eq(k).eig(3));
  v = gsorth([real(eq(k).evec(:,3)) imag(eq(k).evec(:,3)) real(eq(k).evec(:,6))]);
  for delta = [0:0.05:ere*period 0.103 0.1028],
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmedt(L, tend, a0, h, 4); av = [av; v'*aa]; end
  figure(1); set(gcf,'pos',[100 200 700 500]); clf;
  ax1 = axes('pos',[0.1 0.12 0.83 0.82]);
    plot3(av(1:3:end-5,:)',av(2:3:end-4,:)',av(3:3:end-3,:)','-','color',[.5 .5 .5]);
    hold on; grid on; axis equal; view(-110,30); axis off;
    plot3(av(end-2,1:300)',av(end-1,1:300)',av(end,1:300)','k-','linewidth',2);
    e2 = v'*aa(:,end); plot3(e2(1),e2(2),e2(3),'k.','markersize',15);
    e3 = v'*eq(3).a; plot3(e3(1),e3(2),e3(3),'k.','markersize',15);
    text(-0.4, 0.3, 0.52,'E2','fontsize',15);    
    text(-0.4,-0.75, 0.35,'E3','fontsize',15);    
    set(gcf,'paperpos',[3 6 15 11]); wysiwyg;
    print -depsc2 cvitanovic/kse22_E1_UM2_only.eps

%% E1 2nd unstable manifold (with orbits in color)
  clear;  load kse22orbits;  k = 1;  h = 0.1;  tend = 220;  av = [];
  ere = real(eq(k).eig(3));  period = 2*pi/imag(eq(k).eig(3));
  v = gsorth([real(eq(k).evec(:,3)) imag(eq(k).evec(:,3)) real(eq(k).evec(:,6))]);
  for delta = [0:0.06:ere*period 0.8 0.098 0.1028],
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
%    [tt, aa] = ksfmedt(L, tend, a0, h, 2); av = [av; v'*aa];
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2); av = [av; v'*aa];
    if delta == 0.8, aa1 = aa; end
    if delta == 0.098, aa2 = aa; end, end,
  figure(1); set(gcf,'pos',[5 400 420 300],'paperpos',[5 8 12 9]); clf;
  ax1 = axes('pos',[0.06 0.08 0.92 0.90]);
    plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','-','color',[.5 .5 .5]);
    hold on; grid on; axis equal; view(-110,30);
    plot3(av(end-8,1:2:end)',av(end-7,1:2:end)',av(end-6,1:2:end)','.-','color',[0 .8 0]);
    plot3(av(end-5,:)',av(end-4,:)',av(end-3,:)','.-','color','r');
    plot3(av(end-2,:)',av(end-1,:)',av(end,:)','.-','color','b');
    text(-0.4,-0.1,-0.2,'A','fontsize',14,'color',[0 .8 0]);
    text(-0.4,-0.5,-0.3,'B','fontsize',14,'color','r');
    text(-0.4,-0.5, 0.6,'C','fontsize',14,'color','b');
    text(-0.5, 0.3, 0.15,'$E_1$','fontsize',14,'interp','latex');    
    text(-0.4, 0.3, 0.52,'$E_2$','fontsize',14,'interp','latex');    
    text(-0.4,-0.75, 0.35,'$E_3$','fontsize',14,'interp','latex');
    e2 = v'*aa(:,end); plot3(e2(1),e2(2),e2(3),'k.','markersize',15);
    e3 = v'*eq(3).a; plot3(e3(1),e3(2),e3(3),'k.','markersize',15);

    if 0, phase = 0:0.01:1;  ve = eq(2).a(1:2:end)+1i*eq(2).a(2:2:end);
    ve = repmat(ve,1,length(phase)).*exp(-2i*pi*(1:31)'*phase);
    ae = zeros(size(aa,1),length(phase));
    ae(1:2:end,:) = real(ve); ae(2:2:end,:) = imag(ve); ve = v'*ae;
    plot3(ve(1,:),ve(2,:),ve(3,:),'m-','linewidth',1.5);
    phase = 0:0.01:1;  ve = eq(3).a(1:2:end)+1i*eq(3).a(2:2:end);
    ve = repmat(ve,1,length(phase)).*exp(-2i*pi*(1:31)'*phase);
    ae = zeros(size(aa,1),length(phase));
    ae(1:2:end,:) = real(ve); ae(2:2:end,:) = imag(ve); ve = v'*ae;
    plot3(ve(1,:),ve(2,:),ve(3,:),'k-','linewidth',1.5); end

    xlabel('$v_1$','fontsize',14,'interp','latex','pos',[0 1.0 -0.5]);
    ylabel('$v_2$','fontsize',14,'interp','latex','pos',[-0.5 0.0 -0.7]);
    zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex','pos',[1.0 0.85 0]);
    print -depsc2 ks22_E1_plane2_manifold_c.eps
%%
  figure(2); clf; set(gcf,'pos',[450 400 500 300],'paperpos',[5 8 12 8]);
  ax2 = axes('pos',[0.10 0.13 0.23 0.75]);
    a0 = eq(k).a + 1e-4.*exp(0.8).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 15);  [x, uu] = ksfm2real(aa, L); 
    pcolor(x,tt,uu'); caxis([-3 3]); shading flat;
    ylabel('$t$','rotat',0,'fontsize',14,'interp','latex'); 
    xlabel('$x$','fontsize',14,'interp','latex');
    ht = title('Orbit A');  set(ht,'color',[0 .8 0],'fontsize',15);
  ax3 = axes('pos',[0.41 0.13 0.23 0.75]);
    a0 = eq(k).a + 1e-4.*exp(0.098).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 15);  [x, uu] = ksfm2real(aa(1:30,:), L); 
    pcolor(x,tt(1:2:end),uu(:,1:2:end)'); caxis([-3 3]); shading flat;
    xlabel('$x$','fontsize',14,'interp','latex');
    ht = title('Orbit B');  set(ht,'color','r','fontsize',15);
  ax4 = axes('pos',[0.72 0.13 0.23 0.75]);
    a0 = eq(k).a + 1e-4.*exp(0.1028).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 15);  [x, uu] = ksfm2real(aa(1:30,:), L);
    pcolor(x,tt(1:2:end),uu(:,1:2:end)'); caxis([-3 3]); shading flat;
    xlabel('$x$','fontsize',14,'interp','latex');
    ht = title('Orbit C');  set(ht,'color','b','fontsize',15);
    print -depsc ks22_E1_plane2_orbits_c.eps

%% E2 unstable manifold and orbits    
  clear;  load kse22orbits;  k = 2;  h = 0.1;  tend = 150;  av = [];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
  for delta = [0:0.03:1, 0.74 0.066 0.0693]*ere*period,
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 3); av = [av; v'*aa];
    if delta == 0.74*ere*period, aa1 = aa; end
    if delta == 0.066*ere*period, aa2 = aa; end, end
  figure(1); set(gcf,'pos',[5 400 420 300],'paperpos',[5 8 12 9]); clf;
  ax1 = axes('pos',[0.12 0.10 0.83 0.88]);
    plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','-','color',[.5 .5 .5]);
    hold on; grid on; axis equal; view(27,30);
    plot3(av(end-8,:)',av(end-7,:)',av(end-6,:)','.-','color',[0 .8 0]);
    plot3(av(end-5,:)',av(end-4,:)',av(end-3,:)','r.-');
    plot3(av(end-2,:)',av(end-1,:)',av(end,:)','b.-');
    text(0.5, 0.1,-0.5,'A','fontsize',14,'color',[0 .8 0]);
    text(1.0, 0.0,-0.5,'B','fontsize',14,'color','r');
    text(1.0, 0.0, 0.3,'C','fontsize',14,'color','b');
    text(-0.2, 0.3,-0.7,'$E_2$','fontsize',14,'interp','latex');    
%    text(-0.4, 0.3,-0.15,'$\mathbf{g}(L/4)E_2$','fontsize',14,'interp','latex');    
    text(-0.6, 0.3, 0.2,'$\tau_{1/4}E_2$','fontsize',14,'interp','latex');    
    text(1.1, 0.0,-0.1,'$E_3$','fontsize',14,'interp','latex');
    e2 = v'*aa(:,end); plot3(e2(1),e2(2),e2(3),'k.','markersize',16);
    e3 = v'*eq(3).a; plot3(e3(1),e3(2),e3(3),'k.','markersize',16);

    if 0, phase = 0:0.01:1;  ve = eq(3).a(1:2:end)+1i*eq(3).a(2:2:end);
    ve = repmat(ve,1,length(phase)).*exp(-2i*pi*(1:31)'*phase);
    ae = zeros(size(aa,1),length(phase));
    ae(1:2:end,:) = real(ve); ae(2:2:end,:) = imag(ve); ve = v'*ae;
    plot3(ve(1,:),ve(2,:),ve(3,:),'k-','linewidth',1.5); end

    xlabel('$v_1$','fontsize',14,'interp','latex','pos',[0 -.7 -.8]); 
    ylabel('$v_2$','fontsize',14,'interp','latex','pos',[1.3 0 -.8]); 
    zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex','pos',[-1.4 -.7 0]);
    arr1 = annotation(1,'arrow',[0.458 0.5347],[0.7372 0.5791]);
    print -depsc2 ks22_E2_manifold_c.eps
  figure(2); clf; set(gcf,'pos',[450 400 500 300],'paperpos',[5 8 12 8]);
  ax2 = axes('pos',[0.10 0.13 0.24 0.75]);
    a0 = eq(k).a + 1e-4.*exp(0.74*ere*period).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 12);  [x, uu] = ksfm2real(aa(1:30,:), L); 
    pcolor(x,tt,uu'); caxis([-3 3]); shading flat;
    ylabel('$t$','rotat',0,'fontsize',14,'interp','latex'); 
    xlabel('$x$','fontsize',14,'interp','latex');
    ht = title('Orbit A');  set(ht,'color',[0 .8 0],'fontsize',15);
  ax3 = axes('pos',[0.41 0.13 0.24 0.75]);
    a0 = eq(k).a + 1e-4.*exp(0.066*ere*period).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 12);  [x, uu] = ksfm2real(aa(1:30,:), L); 
    pcolor(x,tt,uu'); caxis([-3 3]); shading flat;
    xlabel('$x$','fontsize',14,'interp','latex');
    ht = title('Orbit B');  set(ht,'color','r','fontsize',15);
  ax4 = axes('pos',[0.72 0.13 0.24 0.75]);
    a0 = eq(k).a + 1e-4.*exp(0.0693*ere*period).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 12);  [x, uu] = ksfm2real(aa(1:30,:), L);
    pcolor(x,tt,uu'); caxis([-3 3]); shading flat;
    xlabel('$x$','fontsize',14,'interp','latex');
    ht = title('Orbit C');  set(ht,'color','b','fontsize',15);
    print -depsc2 ks22_E2_orbits_c.eps

%% E3 unstable manifold and orbits    
  clear; load kse22orbits; k = 3; h = 0.1; delta = 1e-4; tend = 110; av = [];
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) eq(k).evec(:,4)]);
  for phi = (0.1:0.1:6)*pi/3,
    a0 = eq(k).a + delta.*(v(:,1).*cos(phi)+v(:,2).*sin(phi));
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);  av = [av; v'*aa]; end
  figure(1); set(gcf,'pos',[5 400 420 430],'paperpos',[5 8 11 11]); clf;
  ax1 = axes('pos',[0.12 0.08 0.83 0.90]);
    plot3(av(1:3:end-2,:)',av(2:3:end-1,:)',av(3:3:end,:)','-','color',[.5 .5 .5]); 
    hold on; grid on;
    text(.5, .5,-1.,'A','fontsize',14,'color',[0 .8 0]);
    text(1., .85, -.3,'B','fontsize',14,'color','r');
    text(1., .25,-.55,'C','fontsize',14,'color','b');
    text(.1, .3, .25,'$E_2$','fontsize',14,'interp','latex');    
    text(1.1, .6,-.4,'$E_3$','fontsize',14,'interp','latex');
    xlabel('$v_1$','fontsize',14,'interp','latex','pos',[0 .8 -1.2]); 
    ylabel('$v_2$','fontsize',14,'interp','latex','pos',[.9 .1 -1.0]); 
    zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex','pos',[.7 -.65 -.5]);
  figure(2); clf; set(gcf,'pos',[450 400 500 300],'paperpos',[5 8 12 8]);
  ax2 = axes('pos',[0.10 0.13 0.24 0.75]);
    phi = atan(v(1,2)./v(1,1)) - pi/2 + 2*pi/3;
    a0 = eq(k).a + delta.*(v(:,1).*cos(phi)+v(:,2).*sin(phi));
    [tt, aa] = ksfmedt(L, 120, a0, h, 10);  [x, uu] = ksfm2real(aa(1:30,:), L);  av = v'*aa; 
    pcolor(x,tt,uu'); caxis([-3 3]); shading flat;
    ylabel('$t$','rotat',0,'fontsize',14,'interp','latex'); 
    xlabel('$x$','fontsize',14,'interp','latex');
    ht = title('Orbit A');  set(ht,'color',[0 .8 0],'fontsize',15);
    axes(ax1); plot3(av(1,:),av(2,:),av(3,:),'.-','color',[0 .8 0]);
  figure(2); ax3 = axes('pos',[0.41 0.13 0.24 0.75]);
    phi = atan(v(1,2)./v(1,1)) + pi/2;
    a0 = eq(k).a + delta.*(v(:,1).*cos(phi)+v(:,2).*sin(phi));
    [tt, aa] = ksfmetd2(a0, L, h, 170, 10);  [x, uu] = ksfm2real(aa(1:30,:), L);  av = v'*aa; 
    pcolor(x,tt,uu'); caxis([-3 3]); shading flat;
    xlabel('$x$','fontsize',14,'interp','latex');
    ht = title('Orbit B');  set(ht,'color','r','fontsize',15);
    axes(ax1); plot3(av(1,:),av(2,:),av(3,:),'r.-');
  figure(2); ax4 = axes('pos',[0.72 0.13 0.24 0.75]);
    phi = atan(v(1,2)./v(1,1)) - pi/2;
    a0 = eq(k).a + delta.*(v(:,1).*cos(phi)+v(:,2).*sin(phi));
    [tt, aa] = ksfmetd2(a0, L, h, 120, 10);  [x, uu] = ksfm2real(aa(1:30,:), L);  av = v'*aa;
    pcolor(x,tt,uu'); caxis([-3 3]); shading flat;
    xlabel('$x$','fontsize',14,'interp','latex');
    ht = title('Orbit C');  set(ht,'color','b','fontsize',15);
    axes(ax1); plot3(av(1,:),av(2,:),av(3,:),'b.-'); 
    if 1, phase = 0:0.01:1;  ve = eq(2).a(1:2:end)+1i*eq(2).a(2:2:end);
      ve = repmat(ve,1,length(phase)).*exp(-1i*pi*(1:31)'*phase);
      ae = zeros(size(aa,1),length(phase));
      ae(1:2:end,:) = real(ve); ae(2:2:end,:) = imag(ve); ve = v'*ae;
      plot3(ve(1,:),ve(2,:),ve(3,:),'k-','linewidth',1.5); end
    axis equal; view(120,30);  set(gca,'zlim',[-1.1 0.1]);
%  print(1,'-depsc2','ks22_E3_manifold.eps');
  print(2,'-depsc2','ks22_E3_orbits_c.eps');

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

%% Plot RPOs in Delta vs T with Delta > 0 and Lyapunov exponents
  clear; load kse22orbits;  tmax = 200;
  trpo = [rpo.T];  phrpo = [rpo.d];
  figure(1); clf; set(gcf,'pos',[200 500 550 200]);
  ipo = find(phrpo < 1e-4); disp(length(ipo));
  trpo(ipo) = [];  phrpo(ipo) = [];
  plot(trpo,phrpo,'k.','markersize',8);
  set(gca,'pos',[0.09 0.20 0.89 0.75],'xtick',0:40:tmax);
  axis([0 tmax 0 L/2]);  xlabel('$T$','fontsize',15,'interp','latex');
  ylabel('$\ell$','fontsize',15,'interp','latex','rotat',0); 
  orient portrait;  set(gcf,'paperpos',[3.0 12.0 15.0 5.5]);
  print -deps ks22_rpos_Tdelta.eps
  
%% Plot largest Lyapunov exponent for RPOs and POs  
  clear; load kse22orbits;  tmax = 200;
  trpo = [rpo.T]; eigrpo = [rpo.eig];
  figure(2); clf; set(gcf,'pos',[200 500 550 200],'paperpos',[3.0 12.0 15.0 5.5]);
  plot(trpo,log(abs(eigrpo(1,:)))./trpo,'k.','markersize',8);
  Lyap = 0.048; hold on; plot([0 tmax],[Lyap Lyap],'r-');
  set(gca,'pos',[0.09 0.20 0.89 0.75],'xtick',0:40:tmax);
  axis([0 tmax 0 0.37]);  
  xlabel('$T$','fontsize',15,'interp','latex');
  ylabel('$s_1$','fontsize',15,'interp','latex','rotat',0,'pos',[-13 0.15]);
%  ylabel('Lyapunov exponent','fontsize',12,'interp','latex'); 
  orient portrait;  set(gcf,'paperpos',[3.0 12.0 15.0 5.5]);
  print -deps ks22_rpos_lyap.eps
  

%% Plot proximity of chaotic orbit to RPOs
  clear; load ks22prox04;   % time average prox:
  iave = 30;  aprox = prox;  cnt = 1;
  for ia = 1:iave,  cnt = cnt + 2; 
    aprox = aprox + prox(:,[ia+1:-1:2 1:end-ia ]); 
    aprox = aprox + prox(:,[ia+1:end end-1:-1:end-ia]); end
  aprox = aprox./cnt;

  [sprox,iprox] = sort(aprox);  kprox = iprox - 7;
  ip = find(mod(kprox,2) == 1);  im = find(mod(kprox,2) == 0);
  kprox(ip) = round((kprox(ip)+1)/2);  kprox(im) = -round(kprox(im)/2);
  
  figure(1); clf; set(gcf,'pos',[5 465 988 480]);
  axes('pos',[0.0525 0.7522 0.9300 0.1889]);
  [x, uu] = ksfm2real(aa, 22, 64);
  pcolor(tt,x,uu); caxis([-3 3]); shading flat;  ylabel('x','rotat',0);

  axes('pos',[0.0525 0.1100 0.9300 0.5800]);
  plot(tt,kprox(1,:),'k.'); set(gca,'ylim',[-72 72]);
  xlabel('Time'); ylabel('RPO index');
  
%% Plotting E2 -> E2 orbits with the unstable manifold of E2
  clear;  load kse22orbits;  irpo = [15 17 19 29 35 45 48 63]; av = [];
  tend = 150; rcol = ['b';'r';'k';'g';'c';'m';'y'];  figure(1); clf;
  ere = real(eq(2).eig(1));  period = 2*pi/imag(eq(2).eig(1));
  v = gsorth([real(eq(2).evec(:,1)) imag(eq(2).evec(:,1)) real(eq(2).evec(:,7))]);
  for delta = (0:0.05:1.05).*ere*period,
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

%% Plotting Relative Equilibria (RE)
  clear;  load kse22orbits; nx = 46; np = 8; tend = 70;
  [x, u1] = ksfm2real(re(1).a, L);
  figure(1); clf; set(gcf,'position',[400 200 400 440],'paperpos',[5 9 10 13]);
  ax1 = axes('pos',[0.10 0.62 0.39 0.34]);  plot(x,u1,'b-','linewidth',1);
  hold on; plot([x([1 end]) repmat(-10:5:10,2,1)],[[0; 0] repmat([-1; 1],1,5)*0.05],'k-');
  set(gca,'xlim',x([1 end]),'xtick',-10:5:10,'ylim',[-2 3]);
%  title('TW$_1$: $c = 0.737$','fontsize',14,'interp','latex');
  text(3,2.5,'TW$_{+1}$','fontsize',14,'interp','latex');
  xlabel('$x$','fontsize',14,'interp','latex'); 
  ylabel('$u$','rotat',0,'fontsize',14,'interp','latex');
%  print -depsc2 ks22_TW1_profile.eps
  figure(2); clf; set(gcf,'position',[400 200 400 440],'paperpos',[5 9 10 13]);
  a0 = re(1).a;  a0(1) = a0(1) + 3e-3;
  [tt, aa] = ksfmetd2(a0, L, h, tend, np);  [x, uu] = ksfm2real(aa(1:nx,:), L);
  ax2 = axes('pos',[0.10 0.10 0.39 0.46]); pcolor(x,tt,uu'); caxis([-3 3]); shading flat;
  xlabel('$x$','fontsize',14,'interp','latex'); 
  ylabel('$t$','rotat',0,'fontsize',12,'interp','latex');
  set(gca,'ytick',0:20:60,'xtick',-10:5:10);
  print -depsc2 ks22_TW1_orbit_c.eps
  figure(3); clf; set(gcf,'position',[400 200 400 440],'paperpos',[5 9 10 13]);
  [x, u2] = ksfm2real(re(2).a, L);
  ax3 = axes('pos',[0.57 0.62 0.39 0.34]);  plot(x,u2,'b-','linewidth',1);  
  hold on; plot([x([1 end]) repmat(-10:5:10,2,1)],[[0; 0] repmat([-1; 1],1,5)*0.05],'k-');
  set(gca,'xlim',x([1 end]),'xtick',-10:5:10,'ylim',[-2 3]);
%  title('TW$_2$: $c = 0.350$','fontsize',14,'interp','latex');
  text(3,2.5,'TW$_{+2}$','fontsize',14,'interp','latex');
  xlabel('$x$','fontsize',14,'interp','latex'); 
  ylabel('$u$','rotat',0,'fontsize',14,'interp','latex');
%  print -depsc2 ks22_TW2_profile.eps
  figure(4); clf; set(gcf,'position',[400 200 400 440],'paperpos',[5 9 10 13]);
  [tt, aa] = ksfmetd2(re(2).a, L, h, tend, np);  [x, uu] = ksfm2real(aa(1:nx,:), L);
  ax4 = axes('pos',[0.57 0.10 0.39 0.46]); pcolor(x,tt,uu'); caxis([-3 3]); shading flat;
  xlabel('$x$','fontsize',14,'interp','latex'); 
  ylabel('$t$','rotat',0,'fontsize',12,'interp','latex');
  set(gca,'ytick',0:20:60,'xtick',-10:5:10);
  print -depsc2 ks22_TW2_orbit_c.eps
%  print -depsc2 ks22_RE1-2.eps

%% Plot individual orbits (one per figure)
  clear;  load kse22orbits; nfm = 32;
%  ost = [1 3 4 5 6]; phs = [0 0 0 0 0]; te = 100; np = 12; % short orbits
%  ost = [15 18 21 37 60]; phs = [0 0 0 0 0]; te = 100; np = 12; % close to cage
  ost = [2 22 27 35 47 57]; phs = [-1.89 6.54 1.58 0.12 -2.2 -3.18]; te = 100; np = 12; % periodic orbits
  for ist = 1:length(ost),
    a0 = rpo(ost(ist)).a;  tend = rpo(ost(ist)).T;  ph = rpo(ost(ist)).d;
%    if refl, a0r = a0;  a0r(1:2:end) = -a0r(1:2:end);  phr = -ph; end
    eig = rpo(ost(ist)).eig(1:8);
    N = length(a0)+2;  ek = exp((2i*pi/L).*ph.*(1:N/2-1)');
    [tt, aa] = ksfmedt(L, tend, a0, h, np);
%    if refl, [tr, ar] = ksfmedt(L, tend, a0r, h, 1); end,
    if 0, % shift orbit in time   
      [maa, ima] = max(sum(aa.^2)');  %disp([ima size(aa,2)]);
      if ima > 1, tt = [tt(ima:end) tt(2:ima)+tend]-tt(ima); aas = aa(:,2:ima);
        vi = (aas(1:2:end,:)+1i*aas(2:2:end,:)).*repmat(ek,1,ima-1);
        aas(1:2:end,:) = real(vi);  aas(2:2:end,:) = imag(vi);
        aa = [aa(:,ima:end) aas]; end, end
    if 1, phase = phs(ist); % shift orbit's phase
      vi = (aa(1:2:end,:) + 1i*aa(2:2:end,:)).*...
           repmat(exp((2i*pi/L).*phase.*(1:N/2-1)'),1,length(tt));
      aa(1:2:end,:) = real(vi);  aa(2:2:end,:) = imag(vi); end
    tti = tt(1:end-1); aai = aa(:,1:end-1); ne = ceil(te./tend); aav = aa;
    for ie = 1:ne-1,
      vi = (aav(1:2:end,:)+1i*aav(2:2:end,:)).*repmat(ek,1,size(aav,2));
      aav(1:2:end,:) = real(vi);  aav(2:2:end,:) = imag(vi);
      aai = [aai aav(:,1:end-1)];  tti = [tti tt(1:end-1)+ie*tend]; end
    [x, uui] = ksfm2real(aai(1:nfm,:), L);

    fig1 = figure('PaperPosition',[8 8 4 10],'PaperOrient','port',...
                  'PaperSize',[20.98 29.68],'Position',[400 300 180 450]);
    pcolor(x,tti,uui'); caxis([-3 3]); shading flat; hold on;
    plot(x([1 end])*ones(1,ne),[tend;tend]*(1:ne),'w-');
    plot(mod([ph;ph]*(1:ne-1),L)-L/2,[(1:ne-1);(2:ne)]*tend,'w-');
%    title([sprintf('$T=%5.1f,~',tend) '\Delta=' sprintf('%6.2f$',ph)],...
%          'interp','latex','fontsize',12); 
    set(gca,'ticklength',0.7*get(gca,'ticklength'),'tickdir','out',...
      'ylim',[0 te],'ytick',20:20:te,'box','off');
    xlabel('$x$','interp','latex','fontsize',15); 
%    if ist == 1, ylabel('$t$','interp','latex','rotat',0,'fontsize',15); end
    print(fig1,'-depsc2',sprintf('ks22rpo%05.1f-%05.2f.eps',tend,ph));
    pause; delete(fig1);
  end

%% Plot individual orbits (one figure)
  clear;  load kse22orbits;
% ost = [1 2 3 4 5];  phs = [0 -1.89 0 0 0]; te = 100; np = 8;  name = 'Short'; % short orbits
% ost = [15 17 20 31 51];  phs = [0 0 0 0 0]; te = 100; np = 10; name = 'Cage'; % close to cage
 ost = [24 55 61 63 125]; phs = [1.58 6.17 0.98 0.12 -1.16]; te = 140; np = 14;name = 'PO'; % periodic orbits
  nst = length(ost);
  fig1 = figure('PaperPosition',[2 8 3*nst 8],'PaperOrient','port',...
                'PaperSize',[20.98 29.68],'Position',[200 500 130*nst 350]);
  hax = subplots(1,nst,[0.03 0.0 0.05 0.04],[0.04 0.01 0.08 0.04]);
  for ist = 1:length(ost),
    a0 = rpo(ost(ist)).a;  tend = rpo(ost(ist)).T;  ph = rpo(ost(ist)).d;
%    if refl, a0r = a0;  a0r(1:2:end) = -a0r(1:2:end);  phr = -ph; end
    eig = rpo(ost(ist)).eig(1:8);
    N = length(a0)+2;  ek = exp((2i*pi/L).*ph.*(1:N/2-1)');
    [tt, aa] = ksfmedt(L, tend, a0, h, np);
%    if refl, [tr, ar] = ksfmedt(L, tend, a0r, h, 1); end,
    if 1, % shift orbit in time
      [maa, ima] = max(sum(aa.^2)');  %disp([ima size(aa,2)]);
      if ima > 1, tt = [tt(ima:end) tt(2:ima)+tend]-tt(ima); aas = aa(:,2:ima);
        vi = (aas(1:2:end,:)+1i*aas(2:2:end,:)).*repmat(ek,1,ima-1);
        aas(1:2:end,:) = real(vi);  aas(2:2:end,:) = imag(vi);
        aa = [aa(:,ima:end) aas]; end, end
    if 1, phase = phs(ist); % shift orbit's phase
      vi = (aa(1:2:end,:) + 1i*aa(2:2:end,:)).*...
           repmat(exp((2i*pi/L).*phase.*(1:N/2-1)'),1,length(tt));
      aa(1:2:end,:) = real(vi);  aa(2:2:end,:) = imag(vi); end
    tti = tt(1:end-1); aai = aa(:,1:end-1); ne = ceil(te./tend); aav = aa;
    for ie = 1:ne-1,
      vi = (aav(1:2:end,:)+1i*aav(2:2:end,:)).*repmat(ek,1,size(aav,2));
      aav(1:2:end,:) = real(vi);  aav(2:2:end,:) = imag(vi);
      aai = [aai aav(:,1:end-1)];  tti = [tti tt(1:end-1)+ie*tend]; end
    [x, uui] = ksfm2real(aai(1:38,:), L);
    
  axes(hax(ist));
    pcolor(x,tti,uui'); caxis([-3 3]); shading flat; hold on;
    plot(x([1 end])*ones(1,ne),[tend;tend]*(1:ne),'w-');
    plot(mod([ph;ph]*(1:ne-1),L)-L/2,[(1:ne-1);(2:ne)]*tend,'w-');
%    title([sprintf('$T=%5.1f,~',tend) '\Delta=' sprintf('%6.2f$',ph)],...
%          'interp','latex','fontsize',12); 
    title(['(' char(ist+96) ')'],'fontsize',14,'interp','latex');
    set(gca,'ticklength',0.7*get(gca,'ticklength'),'tickdir','out',...
      'ylim',[0 te],'ytick',20:20:te,'box','off');
    xlabel('$x$','interp','latex','fontsize',14); 
    if ist == 1, ylabel('$t$','interp','latex','rotat',0,'fontsize',14); end
  end
  print(fig1,'-depsc2',['ks22rpos' name '.eps']);
  pause; delete(fig1);
 
  
%%  Plot g^2 for Ikeda map (to illustrate 'troughness');
 par = [1.0; 0.9; 0.4; 6.0];
 [x0,y0] = bndrgrid([-.5 1.5 -2.5 1],200000);
 [xx,yy] = ikedaxy(par,10,x0,y0);
 g = (xx-x0).^2 + (yy-y0).^2;
 figure(1); clf; pcolor(x0,y0,g); shading flat; hold on;
 x0 = [0.0; 0.1]; nitr = 20000; npre = 1000;
 x = mapitr('ikeda', par, [nitr npre], x0);
 plot(x(1,:),x(2,:),'w.','markersize',2);

%% KS solution with L = n*sqrt(2)*2*pi (n wavelengths) 
  clear;  nw = 10;  L = nw*sqrt(8)*pi;  N = 128;  h = 0.25;
  a0 = zeros(2*N-2,1);  randn('seed',1); a0(1:10) = 0.1*randn(10,1);
  [tt,aa] = ksfmstp(a0, L, h, 800, 4);  
  tt = tt(101:end); aa = aa(:,101:end);
  [xx,uu] = ksfm2real(aa, L);
  fig1 = figure('position',[430  560  540  320],'paperpos',[3 10 14 8]); clf;
  pcolor((xx+L/2)/(sqrt(8)*pi),tt-100,uu'); shading flat;
  set(gca,'pos',[0.10  0.16  0.86  0.78]);
  ylabel('$t$','fontsize',14,'interp','latex','rotat',0);
  xlabel('$x/(2\pi\sqrt{2})$','fontsize',14,'interp','latex');
%  print -depsc ks_largeL.eps
  
%% KS solution with L = n*sqrt(2)*2*pi (with color bar) 
  clear;  nw = 10;  L = nw*sqrt(8)*pi;  N = 128;  h = 0.25;
  a0 = zeros(2*N-2,1);  randn('seed',1); a0(1:10) = 0.1*randn(10,1);
  [tt,aa] = ksfmstp(a0, L, h, 800, 4);  
  tt = tt(101:end); aa = aa(:,101:end);
  [xx,uu] = ksfm2real(aa, L);
  fig1 = figure('position',[430  560  580  320],'paperpos',[3 10 16 8]); clf;
  pcolor((xx+L/2)/(sqrt(8)*pi),tt-100,uu'); shading flat;
  set(gca,'pos',[0.10  0.16  0.86  0.78]);
  ylabel('$t$','fontsize',14,'interp','latex','rotat',0);
  xlabel('$x/(2\pi\sqrt{2})$','fontsize',14,'interp','latex');
  hcb = colorbar;  set(hcb,'pos',[0.90 0.16 0.035 0.78]);
  set(gca,'pos',[0.10 0.16 0.78 0.78]);
  print -depsc ks_largeL_cbar.eps
  
%% Test different fft/ifft conventions in ksfmstp
  clear;  load ks22uqo125a;  L = 50; h = 0.1;
  a0 = [a0; zeros(32,1)];  tend = 3*tend;
  [tt,aa] = ksfmstp(a0, L, h, tend/h, 1);
  [tt,aa2] = ksfmstp2(a0, L, h, tend/h, 1);
  figure(1); clf; plot(tt,aa,'.-',tt,aa2,'o-');
  figure(2); clf; 
  for ii = 1:size(aa,1); plot(tt',[aa(ii,:); aa2(ii,:)]','.-'); title(num2str(ii)); pause;  end