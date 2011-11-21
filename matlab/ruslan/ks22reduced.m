% Visualization of KS L=22 reduced state space.
% ES 2011-11-11

%% E2 unstable manifold and orbits projection on eigenbasis (based on ks22figs.m) 
  clear;  load kse22orbits;  k = 2;  h = 0.1;  tend = 150;  av = [];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
  vInv=mfinvTraj(repmat(eq(k).a,1,3)+v);
  for delta = [0:0.03:1, 0.74 0.066 0.0693]*ere*period,
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av = [av; vInv'*aa];
    if delta == 0.74*ere*period, aa1 = aa; end
    if delta == 0.066*ere*period, aa2 = aa; end, 
  end
  figure(1); clf; %set(gcf,'pos',[5 400 420 300],'paperpos',[5 8 12 9]); 
% %   ax1 = axes('pos',[0.12 0.10 0.83 0.88]);
    plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','-','color',[.5 .5 .5]);
    hold on; grid on; % view(27,30);
    plot3(av(end-8,:)',av(end-7,:)',av(end-6,:)','.-','color',[0 .8 0]);
    plot3(av(end-5,:)',av(end-4,:)',av(end-3,:)','r.-');
    plot3(av(end-2,:)',av(end-1,:)',av(end,:)','b.-');
% %     text(0.5, 0.1,-0.5,'A','fontsize',14,'color',[0 .8 0]);
% %     text(1.0, 0.0,-0.5,'B','fontsize',14,'color','r');
% %     text(1.0, 0.0, 0.3,'C','fontsize',14,'color','b');
%    text(-0.4, 0.3,-0.15,'$\mathbf{g}(L/4)E_2$','fontsize',14,'interp','latex');    
%     text(-0.6, 0.3, 0.2,'$\tau_{1/4}E_2$','fontsize',14,'interp','latex');    
%     text(1.1, 0.0,-0.1,'$E_3$','fontsize',14,'interp','latex');
    e2 = vInv'*mfinv(eq(k).a); plot3(e2(1),e2(2),e2(3),'k.','markersize',30);
    e3 = vInv'*mfinv(eq(3).a); plot3(e3(1),e3(2),e3(3),'k.','markersize',26);
    text(e2(1)-0.1, e2(2)-0.1, e2(3),'$E_2$','fontsize',14,'interp','latex');    
    text(e3(1)-0.1, e3(2)-0.1, e3(3),'$E_3$','fontsize',14,'interp','latex');
% %     if 0, phase = 0:0.01:1;  ve = eq(3).a(1:2:end)+1i*eq(3).a(2:2:end);
% %     ve = repmat(ve,1,length(phase)).*exp(-2i*pi*(1:31)'*phase);
% %     ae = zeros(size(aa,1),length(phase));
% %     ae(1:2:end,:) = real(ve); ae(2:2:end,:) = imag(ve); ve = v'*ae;
% %     plot3(ve(1,:),ve(2,:),ve(3,:),'k-','linewidth',1.5); end    
    xlabel('$v_1$','fontsize',14,'interp','latex');
    ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
    zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);
    xlim([-0.1 1.6]);
    ylim([-0.4 1.2]);
    zlim([-0.02 0.08]);
    view(-4.5,76);
% %     set(gca,'FontSize',14)
%     arr1 = annotation(1,'arrow',[0.458 0.5347],[0.7372 0.5791]);
    print -depsc2 ks22_E2_manifold_inv.eps
    
    
%% E2 unstable manifold and orbits using moving frame (based on ks22figs.m) 
  clear;  load kse22orbits;  k = 2;  h = 0.1;  tend = 150;  av = [];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
  vInv=mfinvTraj(v); % Change mfinvTraj manually.
  for delta = [0:0.03:1, 0.74 0.066 0.0693]*ere*period,
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av = [av; vInv'*aa];
    if delta == 0.74*ere*period, aa1 = aa; end
    if delta == 0.066*ere*period, aa2 = aa; end, 
  end
  figure(1); set(gcf,'pos',[5 400 420 300],'paperpos',[5 8 12 9]); clf;
% %   ax1 = axes('pos',[0.12 0.10 0.83 0.88]);
    plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','-','color',[.5 .5 .5]);
    hold on; grid on; % view(27,30);
    plot3(av(end-8,:)',av(end-7,:)',av(end-6,:)','.-','color',[0 .8 0]);
    plot3(av(end-5,:)',av(end-4,:)',av(end-3,:)','r.-');
    plot3(av(end-2,:)',av(end-1,:)',av(end,:)','b.-');
% %     text(0.5, 0.1,-0.5,'A','fontsize',14,'color',[0 .8 0]);
% %     text(1.0, 0.0,-0.5,'B','fontsize',14,'color','r');
% %     text(1.0, 0.0, 0.3,'C','fontsize',14,'color','b');
%    text(-0.4, 0.3,-0.15,'$\mathbf{g}(L/4)E_2$','fontsize',14,'interp','latex');    
%     text(-0.6, 0.3, 0.2,'$\tau_{1/4}E_2$','fontsize',14,'interp','latex');    
%     text(1.1, 0.0,-0.1,'$E_3$','fontsize',14,'interp','latex');
% %     e2 = vInv'*eq(k).a; plot3(e2(1),e2(2),e2(3),'k.','markersize',30);
% %     e3 = vInv'*eq(3).a; plot3(e3(1),e3(2),e3(3),'k.','markersize',26);
% %     text(e2(1)-0.1, e2(2)-0.1, e2(3),'$E_2$','fontsize',14,'interp','latex');    
% %     text(e3(1)-0.1, e3(2)-0.1, e3(3),'$E_3$','fontsize',14,'interp','latex');
% %     if 0, phase = 0:0.01:1;  ve = eq(3).a(1:2:end)+1i*eq(3).a(2:2:end);
% %     ve = repmat(ve,1,length(phase)).*exp(-2i*pi*(1:31)'*phase);
% %     ae = zeros(size(aa,1),length(phase));
% %     ae(1:2:end,:) = real(ve); ae(2:2:end,:) = imag(ve); ve = v'*ae;
% %     plot3(ve(1,:),ve(2,:),ve(3,:),'k-','linewidth',1.5); end    
    xlabel('$v_1$','fontsize',14,'interp','latex');
    ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
    zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);
% %     xlim([-0.1 1.6]);
% %     ylim([-0.4 1.2]);
% %     zlim([-0.02 0.08]);
% %     view(-4.5,76);
% % %     arr1 = annotation(1,'arrow',[0.458 0.5347],[0.7372 0.5791]);
% %     print -depsc2 ks22_E2_manifold_inv.eps

%% E2 unstable manifold and orbits (based on ks22figs.m) 
  clear;  load kse22orbits;  k = 2;  h = 0.1;  tend = 150;  av = [];
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
%   proj=[1,3,4];
proj=[1,4,5];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  for delta = [0:0.03:1, 0.74 0.066 0.0693]*ere*period,
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av = [av; aa(proj(1),:); aa(proj(2),:); aa(proj(3),:);];
    if delta == 0.74*ere*period, aa1 = aa; end
    if delta == 0.066*ere*period, aa2 = aa; end, 
  end
  figure(1); set(gcf,'pos',[5 400 420 300],'paperpos',[5 8 12 9]); clf;
% %   ax1 = axes('pos',[0.12 0.10 0.83 0.88]);
    plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','-','color',[.5 .5 .5]);
    hold on; grid on; % view(27,30);
    plot3(av(end-8,:)',av(end-7,:)',av(end-6,:)','.-','color',[0 .8 0]);
    plot3(av(end-5,:)',av(end-4,:)',av(end-3,:)','r.-');
    plot3(av(end-2,:)',av(end-1,:)',av(end,:)','b.-');
% %     text(0.5, 0.1,-0.5,'A','fontsize',14,'color',[0 .8 0]);
% %     text(1.0, 0.0,-0.5,'B','fontsize',14,'color','r');
% %     text(1.0, 0.0, 0.3,'C','fontsize',14,'color','b');
%    text(-0.4, 0.3,-0.15,'$\mathbf{g}(L/4)E_2$','fontsize',14,'interp','latex');    
%     text(-0.6, 0.3, 0.2,'$\tau_{1/4}E_2$','fontsize',14,'interp','latex');    
%     text(1.1, 0.0,-0.1,'$E_3$','fontsize',14,'interp','latex');
    e2 = mfinv(eq(k).a); plot3(e2(proj(1)),e2(proj(2)),e2(proj(3)),'k.','markersize',30);
    e3 = mfinv(eq(3).a); plot3(e3(proj(1)),e3(proj(2)),e3(proj(3)),'k.','markersize',30);
%     text(e2(1), e2(2)-0.1, e2(3),'$E_2$','fontsize',14,'interp','latex');    
%     text(e3(1), e3(2)-0.1, e3(3),'$E_3$','fontsize',14,'interp','latex');
% %     if 0, phase = 0:0.01:1;  ve = eq(3).a(1:2:end)+1i*eq(3).a(2:2:end);
% %     ve = repmat(ve,1,length(phase)).*exp(-2i*pi*(1:31)'*phase);
% %     ae = zeros(size(aa,1),length(phase));
% %     ae(1:2:end,:) = real(ve); ae(2:2:end,:) = imag(ve); ve = v'*ae;
% %     plot3(ve(1,:),ve(2,:),ve(3,:),'k-','linewidth',1.5); end    
    xlabel('$v_1$','fontsize',14,'interp','latex');
    ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
    zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);
    xlim([-0.1 1.6]);
    ylim([-0.4 1.2]);
    zlim([-0.02 0.08]);
    view(-4.5,76);
%     arr1 = annotation(1,'arrow',[0.458 0.5347],[0.7372 0.5791]);
    print -depsc2 ks22_E2_manifold_inv.eps


%% Shadowing of RPO(70.35)

clear; load ks22f90h25.mat;
h=0.25; N=16; L=22;
np=1;

refpo=48; % Pick reference cycle

% Load minimum distance from reference cycle data
sfile=['data/ks22ppo' num2str(refpo) '_rpo_cr.dat']; 

crdat=load(sfile);

fig1= figure();

semilogy(crdat(:,1),crdat(:,2),'.');
xlabel('i');
ylabel('d_{min}');

% Integrate reference cycle i.c.
a0=ppo(refpo).a;
[tt0, aa0] = ksfmetd(a0, L, h, ppo(refpo).T, np); 

% Pick shadowing cycle
ipo=6;

a0 = rpo(ipo).a;
[tt, aa] = ksfmetd(a0, L, h, rpo(ipo).T, np); 

% Compute distance of points on cycles
[d, aamf, aa0mf, pos]=minDistanceInvPos(aa,aa0);

globmin=find(d==min(d));

% points of minimal distance (to plot if needed)
if size(aa0mf,2)>=size(aamf,2), 
    pnt0mf = aa0mf(:,pos(globmin));
    pntmf = aamf(:,globmin);
else
    pnt0mf = aa0mf(:,globmin);
    pntmf = aamf(:,pos(globmin));
end

% difference vector at minimal distance
difvmf=pntmf-pnt0mf;

% pick projection for plot
proj=[ 1, 3, 4];

fig2= figure();

plot3(aamf(proj(1),:),aamf(proj(2),:),aamf(proj(3),:),'.-');

hold on;

plot3(aa0mf(proj(1),:),aa0mf(proj(2),:),aa0mf(proj(3),:),'r.-');
% plot3(pnt0mf( proj(1)), pnt0mf(proj(2)), pnt0mf(proj(3)),'ks');
% plot3(pntmf( proj(1)), pntmf(proj(2)), pntmf(proj(3)),'ko');
%%%%

%%% plot alternative shadowing ppo
ipo1=9;
a0=rpo(ipo1).a;
[tt, aa1] = ksfmetd(a0, L, h, rpo(ipo1).T, np); 

[d1, aamf1, aa0mf, pos1]=minDistanceInvPos(aa1,aa0);

plot3(aamf1(proj(1),:),aamf1(proj(2),:),aamf1(proj(3),:),'g-', 'LineWidth', 2);
%%%%
% %%% plot alternative ppo
% ipo2=7;
% a0=ppo(ipo2).a;
% [tt, aa2] = ksfmetd(a0, L, h, ppo(ipo2).T, np); 
% 
% [d, aamf2, aa0mf]=minDistanceInv(aa2,aa0);
% 
% plot3(aamf2(proj(1),:),aamf2(proj(2),:),aamf2(proj(3),:),'k-', 'LineWidth', 2);
% %%%%

xlabel('\beta_1');
ylabel('\beta_2');
zlabel('\gamma_2');

l2=['T_p=' num2str(ppo(refpo).T,4)];
l1=['T_p=' num2str(rpo(ipo).T,4)];
l3=['T_p=' num2str(rpo(ipo1).T,4)];
% l4=['T_p=' num2str(ppo(ipo2).T,4)];
legend(l1,l2,l3);

% use [az, el]=view; 

view(27.5,22)

saveas(fig2, ['ks22ppoT', '7035', 'shad.png']);
%%%%%

% Compute Floquet exponents (from stored multipliers)
fexp_refpo = [log(abs(ppo(refpo).e(1)))/ppo(refpo).T, log(abs(ppo(refpo).e(2)))/ppo(refpo).T];
fexp_ipo = [log(abs(rpo(ipo).e(1)))/rpo(ipo).T, log(abs(rpo(ipo).e(2)))/rpo(ipo).T];
fexp_ipo1 = [log(abs(rpo(ipo1).e(1)))/rpo(ipo1).T, log(abs(rpo(ipo1).e(4)))/rpo(ipo1).T];


%%%%%

%%% Minimum angle along PO(70.35)

% clear; load ks22f90h25angl.mat;
% h=0.25; N=16; L=22;
% np=1;
% refpo=48;

angl=load('ks22ppo_angl_48.dat');

minanglpos=find( angl(:,1) == min(angl(:,1)));

%%%% Plot cycles in phase space, color-code distance from longer orbit.
%%%% Not very clear as it stands
%
% fig3= figure();
% 
% plot3(aa0mf(proj(1),:),aa0mf(proj(2),:),aa0mf(proj(3),:),'.-');
% 
% hold on;
% 
% dmin=min(d);
% dmax=max(d);
% 
% for i=1:size(aamf,2),
%     plot3(aamf(proj(1),i),aamf(proj(2),i),aamf(proj(3),i), '.', 'Color', [d(i)/dmax 0 0]);
% end
% 
% 
% d1min=min(d1);
% d1max=max(d1);
% 
% for i=1:size(aamf1,2),
%     plot3(aamf1(proj(1),i),aamf1(proj(2),i),aamf1(proj(3),i), '.', 'Color', [0 d1(i)/d1max 0]);
% end
% 
% 
% plot3(aa0mf(proj(1),minanglpos),aa0mf(proj(2),minanglpos),aa0mf(proj(3),minanglpos),'ks');

fig4=figure();
hold on;

[AX,H1,H2] =plotyy(pos*h,d,tt0,angl(:,1)/pi);

set(H1,'LineStyle','v', 'Color', [0 0 0.8]);
set(H2,'LineStyle','-', 'Color', [0.8 0 0], 'LineWidth', 2);

plot(pos1*h,d1,'x','Color', [0 0.8 0]);

set(AX(1),'YColor', 'k', 'FontSize', 11); 
set(AX(2),'YColor', 'k', 'FontSize', 11); 
set(get(AX(1),'Ylabel'),'String','d','FontSize', 12);
set(get(AX(2),'Ylabel'),'String','\theta_{1,2}','FontSize', 12); 
xlabel('t');

l1=['d(' num2str(rpo(ipo).T,4) ', ' num2str(ppo(refpo).T,4) ')'];
l3=['\theta_{1,2} for PO(' num2str(ppo(refpo).T,4) ')'];
l2=['d(' num2str(rpo(ipo1).T,4) ', ' num2str(ppo(refpo).T,4) ')'];

saveas(fig4, ['ks22ppoT', '7035', 'angl_dist.png']);
saveas(fig4, ['ks22ppoT', '7035', 'angl_dist.pdf']);
saveas(fig4, ['ks22ppoT', '7035', 'angl_dist.eps']);


