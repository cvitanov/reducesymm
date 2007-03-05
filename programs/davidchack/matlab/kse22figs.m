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
  print -depsc2 cvitanovic/kse22_E1_chaos_color.eps
  
%% E1 2nd unstable manifold
  clear;  load kse22orbits;  k = 1;  h = 0.1;  tend = 220;  av = [];
  ere = real(eq(k).eig(3));  period = 2*pi/imag(eq(k).eig(3));
  v = gsorth([real(eq(k).evec(:,3)) imag(eq(k).evec(:,3)) real(eq(k).evec(:,6))]);
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
%    xlabel('v_1'); ylabel('v_2'); zlabel('v_3','rotat',0);
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
    print -depsc2 cvitanovic/kse22_E1_UM2.eps
    print -dpdf  cvitanovic/kse22_E1_UM2.pdf

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
