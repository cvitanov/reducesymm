% Visualization of KS L=22 reduced state space.
% ES 2011-11-11

%% E2 unstable manifold projection on eigenbasis (based on ks22figs.m) 
  clear;  load kse22orbits;  k = 2;  h = 0.1;  tend = 150;  av = [];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,4))]);
  vInv=mfinvTraj(repmat(eq(k).a,1,3)+v);
  for delta = [0:0.03:1, 0.74 0.066 0.0693]*ere*period,
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
    if delta == 0.74*ere*period, aa01 = aa; end
    if delta == 0.066*ere*period, aa02 = aa; end, 
    if delta == 0.0693*ere*period, aa03 = aa; end, 
%   [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av = [av; vInv'*aa];
  end
  figure(1); clf; %set(gcf,'pos',[5 400 420 300],'paperpos',[5 8 12 9]); 
% %   ax1 = axes('pos',[0.12 0.10 0.83 0.88]);
    plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','-','color',[.5 .5 .5]);
    hold on; grid on; % view(27,30);
    plot3(av(end-8,:)',av(end-7,:)',av(end-6,:)','.-','color',[0 .8 0]);
    plot3(av(end-5,:)',av(end-4,:)',av(end-3,:)','r.-');
    plot3(av(end-2,:)',av(end-1,:)',av(end,:)','b.-');
    e2 = vInv'*mfinvSO2(eq(2).a); plot3(e2(1),e2(2),e2(3),'k.','markersize',30);
    e3 = vInv'*mfinvSO2(eq(3).a); plot3(e3(1),e3(2),e3(3),'k.','markersize',26);
    text(e2(1)-0.04, e2(2)-0.01, e2(3),'$E_2$','fontsize',14,'interp','latex');    
    text(e3(1)+0.06, e3(2)+0.06, e3(3)-0.02,'$E_3$','fontsize',14,'interp','latex');
    xlabel('$v_1$','fontsize',14,'interp','latex');
    ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
    zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);
%     xlim([-0.1 1.6]);
%     ylim([-0.4 1.2]);
%     zlim([-0.02 0.08]);
     view( -159.5, 56);
     axis tight;

     
% %     set(gca,'FontSize',14)
%     arr1 = annotation(1,'arrow',[0.458 0.5347],[0.7372 0.5791]);
    print -depsc2 ks22_E2_manifold_inv.eps
    
    figure(2);clf;
    difs=zeros(size(aa,2),1);
    difs0=zeros(size(aa,2),1);
    for i=1:size(difs,1),
        difs(i)=norm(mfinvSO2(aa03(:,i))-mfinvSO2(eq(k).a));        
        difs0(i)=norm(aa03(:,i)-eq(k).a);        
    end
    semilogy(difs/max(difs),'.'); hold on;
    semilogy(difs0/max(difs0),'r.');
    
    
%% TW1 unstable manifold projection on eigenbasis (based on ks22figs.m) 
  clear; load kse22orbits;  k = 1;  h = 0.1;  tend = 90;  av = [];
  ere = real(re(k).eig(1));  period = 2*pi/imag(re(k).eig(1));
  v = gsorth([real(re(k).evec(:,1)) imag(re(k).evec(:,1)) real(re(k).evec(:,3))]);
  vInv=mfinvTraj(repmat(re(k).a,1,3)+v);
  for delta = [0:0.06:1]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 1);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av = [av; vInv'*aa];
  end
  tend=95; av2=[];
  for delta = [0.495:0.0005:0.5105]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av2 = [av2; vInv'*aa];
  end
  tend=95; av22=[];
  for delta = [0.5105:0.005:0.54]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av22 = [av22; vInv'*aa];
  end

% highlight part of the manifold shadowing shortest rpo.
figure(2);clf; grid on; hold on;
plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','-','color',[.8 .8 .8]);
plot3(av2(1:3:end,:)',av2(2:3:end,:)',av2(3:3:end,:)','-','color',[189/255 183/255 107/255]);
plot3(av22(1:3:end,:)',av22(2:3:end,:)',av22(3:3:end,:)','color', [0 0.8 0]);
r1 = vInv'*mfinvSO2(re(1).a); plot3(r1(1),r1(2),r1(3),'k.','markersize',14);
r2 = vInv'*mfinvSO2(re(2).a); plot3(r2(1),r2(2),r2(3),'k.','markersize',14);
e1 = vInv'*mfinvSO2(eq(1).a); plot3(e1(1),e1(2),e1(3),'k.','markersize',14);
e2 = vInv'*mfinvSO2(eq(2).a); plot3(e2(1),e2(2),e2(3),'k.','markersize',14);
e3 = vInv'*mfinvSO2(eq(3).a); plot3(e3(1),e3(2),e3(3),'k.','markersize',14);
text(r1(1)-0.05, r1(2)-0.05, r1(3),'$TW_1$','fontsize',14, 'interp','latex');    
text(e1(1)-0.05, e1(2)-0.05, e1(3),'$E_1$','fontsize',14,'interp','latex');    
text(e2(1)-0.05, e2(2)-0.05, e2(3),'$E_2$','fontsize',14,'interp','latex');    
% text(e3(1)-0.1, e3(2)-0.1, e3(3),'$E_3$','fontsize',14,'interp','latex');
xlabel('$v_1$','fontsize',14,'interp','latex');
ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);
    
load ks22f90h25t100.mat; np=1; h = 0.1; 
refpo=1;
a0=rpo(refpo).a1;
[tt0, aa0] = ksfmetd(a0, L, h, rpo(refpo).T1, np); 
aa0inv=mfinvTraj(aa0);
av0=vInv'*aa0inv;
plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [218/255 112/255 214/255 ]);



%% Add plot of E2 unstable manifold to previous figure
  h = 0.1;  tend = 100;  av = [];
  ere = real(eq(2).eig(1));  period = 2*pi/imag(eq(2).eig(1));
  v2 = gsorth([real(eq(2).evec(:,1)) imag(eq(2).evec(:,1))  real(eq(2).evec(:,7))]);
  for delta = [0:0.03:1, 0.74 0.066 0.0693]*ere*period,
    a0 = eq(2).a + 1e-4.*exp(delta).*v2(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av = [av; vInv'*aa];
  end
plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','-','color', [	135/255 206/255 250/255]);
plot3(e3(1),e3(2),e3(3),'k.','markersize',14);
text(e3(1)-0.05, e3(2)-0.05, e3(3),'$E_3$','fontsize',14,'interp','latex');
print -depsc2 ks22_TW1_manif_rpo1_E2_manif_inv.eps

%% TW1 unstable manifold projection on eigenbasis 
  clear; load kse22orbits;  k = 1;  h = 0.1;  tend = 90;  av = [];
  ere = real(re(k).eig(1));  period = 2*pi/imag(re(k).eig(1));
  v = gsorth([real(re(k).evec(:,1)) imag(re(k).evec(:,1)) real(re(k).evec(:,3))]);
  vInv=mfinvTraj(repmat(re(k).a,1,3)+v);
    % % Coordinate system using E2, TW+-1 and E3:
    % vInv=gsorth(mfinvTraj([eq(2).a-re(k).a eq(2).a-ksfmRefl(re(k).a) eq(2).a-eq(3).a]));
  for delta = [0:0.06:1]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 1);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av = [av; vInv'*aa];
  end
  
  tend=90; av1=[];
  for delta = [0:0.02:0.10]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av1 = [av1; vInv'*aa];
  end
  tend=90; av2=[];
  for delta = [0.10:0.004:0.158]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av2 = [av2; vInv'*aa];
  end
  tend=95; av3=[];
  for delta = [0.158:0.002:0.18]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av3 = [av3; vInv'*aa];
  end
  tend=100; av4=[];
  for delta = [0.18:0.002:0.21]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av4 = [av4; vInv'*aa];
  end
% visit E2  
  tend=100; av5=[];
  for delta = [0.21:0.01:0.460]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av5 = [av5; vInv'*aa];
  end
% visit E2  
  tend=110; av51=[];
  for delta = [0.3:0.005:0.44]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av51 = [av51; vInv'*aa];
  end

% visit E2 (not so close, contain rpo 5) 
  tend=135; av511=[];
  for delta = [0.441:0.0005:0.446]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av511 = [av511; vInv'*aa];
  end  
  % Reflected TW1 - visit E2  
  tend=110; av51r=[];
  for delta = [0.3:0.005:0.44]*ere*period,
    a0 = ksfmRefl(re(k).a + 1e-4.*exp(delta).*v(:,2));
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av51r = [av51r; vInv'*aa];
  end
% visit E2
  tend=100; av6=[];
  for delta = [0.460:0.001:0.489]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av6 = [av6; vInv'*aa];
  end
  tend=95; av7=[];
  for delta = [0.490:0.001:0.495, 0.485]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av7 = [av7; vInv'*aa];
  end
% shadowing shortest rpo
  tend=107; av81=[]; 
  for delta = [0.5091:0.00001:0.5092]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av81 = [av81; vInv'*aa];
  end
% shadowing shortest rpo
  tend=98; av82=[];
  for delta = [0.509:0.0005:0.51]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av82 = [av82; vInv'*aa];
  end  
  tend=95; av8=[];
  for delta = [0.495:0.0005:0.5105]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av8 = [av8; vInv'*aa];
  end
% 
  tend=95; av9=[];
  for delta = [0.5105:0.005:0.54]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av9 = [av9; vInv'*aa];
  end
  tend=85; av10=[];
  for delta = [0.54:0.04:0.98]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av10 = [av10; vInv'*aa];
  end
  
%   tend=201; av111=[];
%   for delta = 0.4*ere*period,
%     a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
%     [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
% %    [aa, ss] = ksfm_so2(aa, L);
%     aa = mfinvTraj(aa);
%     av111 = [av111; vInv'*aa];
%   end
  
  
figure(1); clf; 
plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','-','color',[.5 .5 .5]);
hold on; grid on; % view(27,30);
plot3(av(end-2,:)',av(end-1,:)',av(end,:)','color',[0 .8 0]);
r1 = vInv'*mfinvSO2(re(1).a); plot3(r1(1),r1(2),r1(3),'k.','markersize',14);
r2 = vInv'*mfinvSO2(re(2).a); plot3(r2(1),r2(2),r2(3),'k.','markersize',14);
e1 = vInv'*mfinvSO2(eq(1).a); plot3(e1(1),e1(2),e1(3),'k.','markersize',14);
e2 = vInv'*mfinvSO2(eq(2).a); plot3(e2(1),e2(2),e2(3),'k.','markersize',14);
e3 = vInv'*mfinvSO2(eq(3).a); plot3(e3(1),e3(2),e3(3),'k.','markersize',14);
text(r1(1)-0.1, r1(2)-0.1, r1(3),'$TW_1$','fontsize',14, 'interp','latex');    
text(e1(1)-0.1, e1(2)-0.1, e1(3),'$E_1$','fontsize',14,'interp','latex');    
text(e2(1)-0.1, e2(2)-0.1, e2(3),'$E_2$','fontsize',14,'interp','latex');    
text(e3(1)-0.1, e3(2)-0.1, e3(3),'$E_3$','fontsize',14,'interp','latex');
xlabel('$v_1$','fontsize',14,'interp','latex');
ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);
set(gca,'FontSize',14)
% arr1 = annotation(1,'arrow',[0.458 0.5347],[0.7372 0.5791]);
% print -depsc2 ks22_TW1_manifold_inv.eps

figure(3);clf; grid off; axis on; hold on;
plot3(av1(1:3:end,:)',av1(2:3:end,:)',av1(3:3:end,:)','-','color',[189/255 183/255 107/255]);
plot3(av2(1:3:end,:)',av2(2:3:end,:)',av2(3:3:end,:)','k-');
plot3(av3(1:3:end,:)',av3(2:3:end,:)',av3(3:3:end,:)','color',[139/255 69/255 19/255]); 
plot3(av4(1:3:end,:)',av4(2:3:end,:)',av4(3:3:end,:)','-','color',[139/255 69/255 19/255]); 
plot3(av5(1:3:end,:)',av5(2:3:end,:)',av5(3:3:end,:)','r-'); % approach E2
plot3(av51(1:3:end,:)',av51(2:3:end,:)',av51(3:3:end,:)','r-'); % approach E2
plot3(av6(1:3:end,:)',av6(2:3:end,:)',av6(3:3:end,:)','r-'); % approach E2
plot3(av7(1:3:end,:)',av7(2:3:end,:)',av7(3:3:end,:)','m-');
plot3(av8(1:3:end,:)',av8(2:3:end,:)',av8(3:3:end,:)','b-'); 
plot3(av9(1:3:end,:)',av9(2:3:end,:)',av9(3:3:end,:)','-', 'color', [0 0.8 0]); 
plot3(av10(1:3:end,:)',av10(2:3:end,:)',av10(3:3:end,:)','-', 'color', [189/255 183/255 107/255]);

plot3(r1(1),r1(2),r1(3),'k.','markersize',20);
plot3(e1(1),e1(2),e1(3),'k.','markersize',20);
plot3(e2(1),e2(2),e2(3),'k.','markersize',20);
plot3(e3(1),e3(2),e3(3),'k.','markersize',20);
text(r1(1)-0.01, r1(2)-0.03, r1(3),'$TW_1$','fontsize',14, 'interp','latex');    
text(e1(1)-0.01, e1(2)-0.03, e1(3),'$E_1$','fontsize',14,'interp','latex');    
text(e2(1)-0.01, e2(2)-0.03, e2(3),'$E_2$','fontsize',14,'interp','latex');    
text(e3(1)-0.01, e3(2)-0.03, e3(3),'$E_3$','fontsize',14,'interp','latex');
xlabel('$v_1$','fontsize',14,'interp','latex');
ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);
set(gca,'FontSize',14)
print -depsc2 ks22_TW1_manifold_inv.eps



figure(4);clf; grid off; axis on; hold on; view(1,30);
plot3(av82(1:3:end,:)',av82(2:3:end,:)',av82(3:3:end,:)','b-'); % shadow rpo 1print -depsc2 ks22_TW1_E2_manif_inv.eps

plot3(av81(1:3:end,:)',av81(2:3:end,:)',av81(3:3:end,:)','-', 'color', [0 0.8 0]); % shadow rpo 1
plot3(r1(1),r1(2),r1(3),'k.','markersize',20);
% plot3(e1(1),e1(2),e1(3),'k.','markersize',20);
% plot3(e2(1),e2(2),e2(3),'k.','markersize',20);
text(r1(1)-0.01, r1(2)-0.01, r1(3)-0.05,'$TW_1$','fontsize',14, 'interp','latex');    
% text(e1(1)-0.01, e1(2)-0.01, e1(3),'$E_1$','fontsize',14,'interp','latex');    
% text(e2(1)-0.01, e2(2)-0.01, e2(3),'$E_2$','fontsize',14,'interp','latex');    
xlabel('$v_1$','fontsize',14,'interp','latex');
ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);
set(gca,'FontSize',14)
load ks22f90h25.mat; np=1; h = 0.1; 
refpo=1;
a0=rpo(refpo).a1;
[tt0, aa0] = ksfmetd(a0, L, h, rpo(refpo).T1, np); 
aa0inv=mfinvTraj(aa0);
av0=vInv'*aa0inv;
plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [218/255 112/255 214/255 ]);
print -depsc2 ks22_TW1_manif_rpo1_inv.eps

figure(5);clf; grid off; axis on; hold on;
% plot3(av5(1:3:end,:)',av5(2:3:end,:)',av5(3:3:end,:)','r-');  % approach E2
plot3(av51(1:3:end,:)',av51(2:3:end,:)',av51(3:3:end,:)','r-'); % approach E2
% % visit E2  (pin heteroclinic)
% ere2 = real(re(k).eig(3));  period2 = 2*pi/imag(re(k).eig(3));
% tend=150; av51h=[]; a51hD=[];
%    for delta2 = (0:0.1:1)*ere2*period2,
%      for delta = [0.383:0.00005:0.384]*ere*period,
%         a0 = re(k).a + 1e-4.*exp(delta).*v(:,2)-1e-4.*exp(delta2).*real(re(k).evec(:,3)); % It was: 
%         [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
%     %    [aa, ss] = ksfm_so2(aa, L);
%         aa = mfinvTraj(aa);
%         dist=zeros(size(aa,2),1);
%         for i=1:size(aa,2),
%             dist(i)=norm(aa(:,i)-mfinvSO2(eq(2).a));    
%         end
%         a51hD = [a51hD, min(dist)];
%         av51h = [av51h; vInv'*aa];
%     end
%    end
% %   for delta = [0.385:0.001:0.390]*ere*period,
% %     a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
% %     [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
% % %    [aa, ss] = ksfm_so2(aa, L);
% %     aa = mfinvTraj(aa);
% %     av51h = [av51h; vInv'*aa];
% %   end
% plot3(av51h(1:3:end,:)',av51h(2:3:end,:)',av51h(3:3:end,:)','g-'); % approach E2
plot3(av511(1:3:end,:)',av511(2:3:end,:)',av511(3:3:end,:)','k-'); % approach E2
plot3(av51r(1:3:end,:)',av51r(2:3:end,:)',av51r(3:3:end,:)','g-'); % approach E1
% plot3(av6(1:3:end,:)',av6(2:3:end,:)',av6(3:3:end,:)','r-'); % approach E1
tend = 100;  avE2 = [];
ere = real(eq(2).eig(1));  period = 2*pi/imag(eq(2).eig(1));
v2 = gsorth([real(eq(2).evec(:,1)) imag(eq(2).evec(:,1))  real(eq(2).evec(:,7))]);
for delta = [0:0.01:1]*ere*period,
    a0 = eq(2).a + 1e-4.*exp(delta).*v2(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 1);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    avE2 = [avE2; vInv'*aa];
  end
plot3(avE2(1:3:end-8,:)',avE2(2:3:end-7,:)',avE2(3:3:end-6,:)','-','color', [	135/255 206/255 250/255]);
r1 = vInv'*mfinvSO2(re(1).a); plot3(r1(1),r1(2),r1(3),'k.','markersize',14);
% r2 = vInv'*mfinvSO2(re(2).a); plot3(r2(1),r2(2),r2(3),'k.','markersize',14);
% e1 = vInv'*mfinvSO2(eq(1).a); plot3(e1(1),e1(2),e1(3),'k.','markersize',14);
e2 = vInv'*mfinvSO2(eq(2).a); plot3(e2(1),e2(2),e2(3),'k.','markersize',22);
e3 = vInv'*mfinvSO2(eq(3).a); plot3(e3(1),e3(2),e3(3),'k.','markersize',14);
text(r1(1)-0.01, r1(2)-0.03, r1(3)-0.05,'$TW_1$','fontsize',14, 'interp','latex');    
% text(e1(1)-0.1, e1(2)-0.1, e1(3),'$E_1$','fontsize',14,'interp','latex');    
text(e2(1)-0.02, e2(2)-0.02, e2(3)+0.05,'$E_2$','fontsize',14,'interp','latex');    
text(e3(1)+0.01, e3(2)+0.02, e3(3),'$E_3$','fontsize',14,'interp','latex');
xlabel('$v_1$','fontsize',14,'interp','latex');
ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);
set(gca,'FontSize',14);
view(13,-4);
axis tight;
print -depsc2 ks22_TW1_E2_manif_inv.eps


figure(6);clf; grid off; axis on; hold on;
plot3(av1(1:3:end,:)',av1(2:3:end,:)',av1(3:3:end,:)','-','color', [0.5 0.5 0.5] );
plot3(av2(1:3:end,:)',av2(2:3:end,:)',av2(3:3:end,:)','color', [0.5 0.5 0.5]);
plot3(av3(1:3:end,:)',av3(2:3:end,:)',av3(3:3:end,:)','color',[0.5 0.5 0.5]); 
plot3(av4(1:3:end,:)',av4(2:3:end,:)',av4(3:3:end,:)','-','color',[0.5 0.5 0.5]); 
plot3(av5(1:3:end,:)',av5(2:3:end,:)',av5(3:3:end,:)','r-'); % approach E1
plot3(av51(1:3:end,:)',av51(2:3:end,:)',av51(3:3:end,:)','r-'); % approach E1
plot3(av6(1:3:end,:)',av6(2:3:end,:)',av6(3:3:end,:)','r-'); % approach E1
plot3(av7(1:3:end,:)',av7(2:3:end,:)',av7(3:3:end,:)', '-', 'color', [0.5 0.5 0.5]);
plot3(av81(1:3:end,:)',av81(2:3:end,:)',av81(3:3:end,:)','b-'); % shadow rpo 1
plot3(av82(1:3:end,:)',av82(2:3:end,:)',av82(3:3:end,:)','-', 'color', [0 0.8 0]); % shadow rpo 1
plot3(av8(1:3:end,:)',av8(2:3:end,:)',av8(3:3:end,:)','-', 'color', [0.5 0.5 0.5]);
plot3(av10(1:3:end,:)',av10(2:3:end,:)',av10(3:3:end,:)','-', 'color', [0.5 0.5 0.5]);
print -depsc2 ks22_TW1_manif_art_inv.eps
plot3(r1(1),r1(2),r1(3),'k.','markersize',20);
% plot3(e1(1),e1(2),e1(3),'k.','markersize',20);
plot3(e2(1),e2(2),e2(3),'k.','markersize',20);
text(r1(1)-0.01, r1(2)-0.03, r1(3)-0.05,'$TW_1$','fontsize',14, 'interp','latex');    
% text(e1(1)-0.1, e1(2)-0.1, e1(3),'$E_1$','fontsize',14,'interp','latex');    
text(e2(1)-0.02, e2(2)-0.02, e2(3),'$E_2$','fontsize',14,'interp','latex');  
xlabel('$v_1$','fontsize',14,'interp','latex');
ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);
set(gca,'FontSize',14);
print -depsc2 ks22_TW1_manif_all_inv.eps

figure(7);clf; grid off; axis on; hold on;
plot3(av1(1:3:end,:)',av1(2:3:end,:)',av1(3:3:end,:)','-','color', [0.5 0.5 0.5] );
plot3(av2(1:3:end,:)',av2(2:3:end,:)',av2(3:3:end,:)','color', [0.5 0.5 0.5]);
plot3(av3(1:3:end,:)',av3(2:3:end,:)',av3(3:3:end,:)','color',[0.5 0.5 0.5]); 
plot3(av4(1:3:end,:)',av4(2:3:end,:)',av4(3:3:end,:)','-','color',[0.5 0.5 0.5]); 
plot3(av5(1:3:end,:)',av5(2:3:end,:)',av5(3:3:end,:)','-','color',[0.5 0.5 0.5]); % approach E1
plot3(av51(1:3:end,:)',av51(2:3:end,:)',av51(3:3:end,:)','-','color',[0.5 0.5 0.5]); % approach E1
plot3(av6(1:3:end,:)',av6(2:3:end,:)',av6(3:3:end,:)','-','color',[0.5 0.5 0.5]); % approach E1
plot3(av7(1:3:end,:)',av7(2:3:end,:)',av7(3:3:end,:)', '-', 'color', [0.5 0.5 0.5]);
plot3(av8(1:3:end,:)',av8(2:3:end,:)',av8(3:3:end,:)','-','color',[0.5 0.5 0.5]);
plot3(av9(1:3:end,:)',av9(2:3:end,:)',av9(3:3:end,:)','-','color',[0.5 0.5 0.5]);
plot3(av10(1:3:end,:)',av10(2:3:end,:)',av10(3:3:end,:)','-','color', [0.5 0.5 0.5]);
plot3(r1(1),r1(2),r1(3),'k.','markersize',20);
% plot3(e1(1),e1(2),e1(3),'k.','markersize',20);
plot3(e2(1),e2(2),e2(3),'k.','markersize',20); 
text(r1(1)-0.01, r1(2)-0.03, r1(3)-0.05,'$TW_1$','fontsize',14, 'interp','latex');    
% text(e1(1)-0.1, e1(2)-0.1, e1(3),'$E_1$','fontsize',14,'interp','latex');    
text(e2(1)-0.02, e2(2)-0.02, e2(3),'$E_2$','fontsize',14,'interp','latex');    
xlabel('$v_1$','fontsize',14,'interp','latex');
ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);
set(gca,'FontSize',14);
% print -depsc2 ks22_TW1_manif_grey_inv.eps


% Add "natural measure" to previous plot
% load ks22f90h25.mat; h = 0.1; 
% % figure(8);
% refpo=62;
% a0=rpo(refpo).a1;
% [tt, aa] = ksfmetd2(a0, L, h, 28*rpo(refpo).T1, 5); 
% aa = mfinvTraj(aa);
% avNM = vInv'*aa;
% plot3(avNM(1,:)',avNM(2,:)',avNM(3,:)','r-');


%% Add rpos and unstable manifold of E2 to previous plot.

load ks22f90h25.mat; np=2; h = 0.1; 
clrs=[ [135/255 100/255 250/255]; [189/255 183/255 107/255]; [139/255 69/255 19/255]; [255/255 165/255 0/255]; [0.1 0.9 0.1] ];
for refpo=1:5,
    a0=rpo(refpo).a1;
    [tt0, aa0] = ksfmetd(a0, L, h, rpo(refpo).T1, np); 
    aa0inv=mfinvTraj(aa0);
    av0=vInv'*aa0inv;
    plot3(av0(1,:)',av0(2,:)',av0(3,:)','-', 'color', clrs(refpo,:), 'LineWidth', 2 ); %[rand(2,1); 0]
end

% figure(); hold on;

tend = 100;  avE2 = [];
ere = real(eq(2).eig(1));  period = 2*pi/imag(eq(2).eig(1));
v2 = gsorth([real(eq(2).evec(:,1)) imag(eq(2).evec(:,1))  real(eq(2).evec(:,7))]);
for delta = [0:0.01:1, 0.74 0.066 0.0693]*ere*period,
    a0 = eq(2).a + 1e-4.*exp(delta).*v2(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 1);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    avE2 = [avE2; vInv'*aa];
  end
plot3(avE2(1:3:end-8,:)',avE2(2:3:end-7,:)',avE2(3:3:end-6,:)','-','color', [	135/255 206/255 250/255]);
e3 = vInv'*mfinvSO2(eq(3).a); plot3(e3(1),e3(2),e3(3),'k.','markersize',14);
text(e3(1)+0.03, e3(2)+0.05, e3(3),'$E_3$','fontsize',14,'interp','latex');
set(gca,'FontSize',14);
axis tight;
print -depsc2 ks22_TW1_E2_manif_rpos_inv.eps

refpo=62;
a0=rpo(refpo).a1;
[tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
aa0inv=mfinvTraj(aa0);
av0=vInv'*aa0inv;
plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [139/255 69/255 19/255]);
refpo=17;
a0=ksfmRefl(rpo(refpo).a1);
[tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
aa0inv=mfinvTraj(aa0);
av0=vInv'*aa0inv;
plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [189/255 183/255 107/255]);
refpo=1;
a0=rpo(refpo).a1;
[tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
aa0inv=mfinvTraj(aa0);
av0=vInv'*aa0inv;
plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [255/255 165/255 0/255]);


%% Plot homoclinic rpo
clear; load kse22orbits;  k = 2;  h = 0.1;
v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
slc=1;    
vInv=mfinvTraj([eq(k).a eq(k).a eq(k).a]+v);
figure(4); clf; view(0,90);
hold on; grid on; % view(27,30);
% r1 = vInv'*mfinvTraj(re(1).a); plot3(r1(1),r1(2),r1(3),'k.','markersize',14);
% e1 = vInv'*mfinvTraj(eq(1).a); plot3(e1(1),e1(2),e1(3),'k.','markersize',14);
e2 = vInv'*mfinvTraj(eq(2).a); plot3(e2(1),e2(2),e2(3),'k.','markersize',24);
e3 = vInv'*mfinvTraj(eq(3).a); plot3(e3(1),e3(2),e3(3),'k.','markersize',20);
% text(r1(1)-0.05, r1(2)-0.05, r1(3),'$TW_1$','fontsize',14,  'interp','latex');    
% text(e1(1)-0.05, e1(2)-0.05, e1(3),'$E_1$','fontsize',14,'interp','latex');    
text(e2(1)-0.05, e2(2)-0.03, e2(3),'$E_2$','fontsize',14,'interp','latex');    
text(e3(1)-0.05, e3(2)+0.03, e3(3),'$E_3$','fontsize',14,'interp','latex');
 xlabel('$v_1$','fontsize',14,'interp','latex');
ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);

% figure(3);clf; grid off; axis on; hold on;

load ks22f90h25.mat; np=1;h = 0.1; 

refpo=62;
a0=rpo(refpo).a1;
[tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
aa0inv=mfinvTraj(aa0);
av0=vInv'*aa0inv;
plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [255/255 165/255 0/255]);

refpo=5;
a0=rpo(refpo).a1;
[tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
aa0inv=mfinvTraj(aa0);
av0=vInv'*aa0inv;
plot3(av0(1,:)',av0(2,:)',av0(3,:)','r.-');


%% Add plot of E2 unstable manifold to previous figure
  k = 2;  h = 0.1;  tend = 100;  av = [];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v2 = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1))  real(eq(k).evec(:,7))]);
  for delta = [0:0.03:1, 0.74 0.066 0.0693]*ere*period,
    a0 = eq(k).a + 1e-4.*exp(delta).*v2(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 1);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av = [av; vInv'*aa];
  end
plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','-','color', [	135/255 206/255 250/255]);
 
print -depsc2 ks22_E2_manif_homo_rpo_inv.eps


%% TW1 unstable manifold (thick) projection on eigenbasis (based on ks22figs.m) 
  clear; load kse22orbits;  k = 1;  h = 0.1;  tend = 90;  av =[]; av2 = []; av3=[];
  ere = real(re(k).eig(1));  period = 2*pi/imag(re(k).eig(1));
  v = gsorth([real(re(k).evec(:,1)) imag(re(k).evec(:,1)) real(re(k).evec(:,3))]);
  v2 = gsorth([real(re(k).evec(:,3)) imag(re(k).evec(:,3)) real(re(k).evec(:,1))]);
  ere2 = real(re(k).eig(3));  period2 = 2*pi/imag(re(k).eig(3));
  vInv=mfinvTraj(repmat(re(k).a,1,3)+v);
   for delta = (0:0.03:1)*ere*period,
        a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
        [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
    %    [aa, ss] = ksfm_so2(aa, L);
        aa = mfinvTraj(aa);
        av = [av; vInv'*aa];
  end
   for delta2 = (0:0.1:1)*ere2*period2,
    for delta = (0:0.03:1)*ere*period,
        a0 = re(k).a + 1e-4.*exp(delta).*v(:,2)+ 1e-4.*exp(delta2).*real(re(k).evec(:,3));
        [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
    %    [aa, ss] = ksfm_so2(aa, L);
        aa = mfinvTraj(aa);
        av2 = [av2; vInv'*aa];
    end
   end
%   for delta2 = [0:0.03:1, 0.06 0.12 0.18]*ere2*period2,
%     a0 = re(k).a + 1e-4.*exp(delta).*v(:,2)+ 1e-4.*exp(delta2).*v(:,3);
%     [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
% %    [aa, ss] = ksfm_so2(aa, L);
%     aa = mfinvTraj(aa);
%     av2 = [av2; vInv'*aa];
%   end
%   for delta2 = [0:0.03:1, 0.06 0.12 0.18]*ere2*period2,
%     a0 = re(k).a + 1e-4.*exp(delta/2).*v(:,2)+ 1e-4.*exp(delta2).*v(:,3);
%     [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
% %    [aa, ss] = ksfm_so2(aa, L);
%     aa = mfinvTraj(aa);
%     av3 = [av3; vInv'*aa];
%   end
  figure(1); clf; %set(gcf,'pos',[5 400 420 300],'paperpos',[5 8 12 9]); 
% %   ax1 = axes('pos',[0.12 0.10 0.83 0.88]);
    plot3(av(1:3:end,:)',av(2:3:end,:)',av(3:3:end,:)','-','color',[.5 .5 .5]);
    hold on; grid on; % view(27,30);
     plot3(av2(1:3:end-8,:)',av2(2:3:end-7,:)',av2(3:3:end-6,:)','r-');
%     plot3(av3(1:3:end-8,:)',av3(2:3:end-7,:)',av3(3:3:end-6,:)','b-');
    r1 = vInv'*mfinvTraj(re(1).a); plot3(r1(1),r1(2),r1(3),'k.','markersize',14);
    e1 = vInv'*mfinvTraj(eq(1).a); plot3(e1(1),e1(2),e1(3),'k.','markersize',14);
    e2 = vInv'*mfinvTraj(eq(2).a); plot3(e2(1),e2(2),e2(3),'k.','markersize',14);
    e3 = vInv'*mfinvTraj(eq(3).a); plot3(e3(1),e3(2),e3(3),'k.','markersize',14);
    text(r1(1)-0.1, r1(2)-0.1, r1(3),'$TW_1$','fontsize',14, 'color', [1 1 1], 'interp','latex');    
    text(e1(1)-0.1, e1(2)-0.1, e1(3),'$E_1$','fontsize',14,'interp','latex');    
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
%     xlim([-0.1 1.6]);
%     ylim([-0.4 1.2]);
%     zlim([-0.02 0.08]);
     view(23.5,10);
% %     set(gca,'FontSize',14)
%     arr1 = annotation(1,'arrow',[0.458 0.5347],[0.7372 0.5791]);
    print -depsc2 ks22_TW1_manifold_thick_inv.eps
    
%% Add rpos to previous plot.

load ks22f90h25.mat; np=1; h = 0.1; 
refpo=1;
a0=rpo(refpo).a1;
[tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
aa0inv=mfinvTraj(aa0);
av0=vInv'*aa0inv;
plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [218/255 112/255 214/255 ]);
% refpo=3;
% a0=rpo(refpo).a1;
% [tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
% aa0inv=mfinvTraj(aa0);
% av0=vInv'*aa0inv;
% plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [189/255 183/255 107/255]);
% refpo=4;
% a0=rpo(refpo).a1;
% [tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
% aa0inv=mfinvTraj(aa0);
% av0=vInv'*aa0inv;
% plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [139/255 69/255 19/255]);
% refpo=6;
% a0=rpo(refpo).a1;
% [tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
% aa0inv=mfinvTraj(aa0);
% av0=vInv'*aa0inv;
% plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [255/255 165/255 0/255]);

print -depsc2 ks22_TW1_manifold_thick_rpo1_inv.eps



%% E2 unstable manifold using slice c_1=0, b_1>0 (based on ks22figs.m) 
  clear;  load kse22orbits;  k = 2;  h = 0.1;  tend = 150;  av = [];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
  vInv=mf1Traj(v);
  for delta = [0:0.03:1, 0.74 0.066 0.0693]*ere*period,
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mf1Traj(aa);
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


%% E2 unstable manifold using slice c_i=0, b_i>0 (based on ks22figs.m) 
  clear;  load kse22orbits;  k = 2;  h = 0.1;  tend = 150;  av = [];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
  % Pick slice:
  slc=1;
  vInv=mfTraj(v,slc);
  for delta = [0:0.03:1, 0.74 0.066 0.0693]*ere*period,
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
%    [aa, ss] = ksfm_so2(aa, L);
     [aa, ss] = mfTraj(aa,slc);
      nj = 100;  tj = zeros(1,nj);
      ij = find(abs(diff(ss))>0.9*pi/slc);
      tj(1,1:length(ij)) = (tt(ij)+tt(ij+1))/2;
      tas = [tj tt];  tas = [tas; [NaN(size(aa,1),nj) aa]]';
      tas = sortrows(tas,1)';  aa = tas(2:end,:);
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

%% TW1 unstable manifold projection on eigenbasis using slice c_i-0, b_i>0 (based on ks22figs.m) 
  clear; load kse22orbits;  k = 1;  h = 0.1;  tend = 90;  av = [];
  ere = real(re(k).eig(1));  period = 2*pi/imag(re(k).eig(1));
  v = gsorth([real(re(k).evec(:,1)) imag(re(k).evec(:,1)) real(re(k).evec(:,3))]);
  slc=1;    
  vInv=mfTraj([re(k).a re(k).a re(k).a]+v,slc);
  for delta = [0:0.03:1, 0.06 0.51 0.54]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
%    [aa, ss] = ksfm_so2(aa, L);
    [aa, ss] = mfTraj(aa,slc);
    nj = 100;  tj = zeros(1,nj);
    ij = find(abs(diff(ss))>0.9*pi/slc);
    tj(1,1:length(ij)) = (tt(ij)+tt(ij+1))/2;
    tas = [tj tt];  tas = [tas; [NaN(size(aa,1),nj) aa]]';
    tas = sortrows(tas,1)';  aa = tas(2:end,:);
    av = [av; vInv'*aa];
  end
  figure(2); clf; %set(gcf,'pos',[5 400 420 300],'paperpos',[5 8 12 9]); 
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
    r1 = vInv'*mfTraj(re(1).a,slc); plot3(r1(1),r1(2),r1(3),'k.','markersize',14);
    [e1 s] = mfTraj(eq(1).a,slc);
    e1 = vInv'*e1; plot3(e1(1),e1(2),e1(3),'k.','markersize',14);
    e2 = vInv'*fmRot(eq(2).a,s); plot3(e2(1),e2(2),e2(3),'k.','markersize',14);
    e2s = vInv'*fmRot(eq(2).a,s+pi/2); plot3(e2(1),e2(2),e2(3),'k.','markersize',14);
    e3 = vInv'*fmRot(eq(3).a,s); plot3(e3(1),e3(2),e3(3),'k.','markersize',14);
    text(r1(1)-0.1, r1(2)-0.1, r1(3),'$TW_1$','fontsize',14, 'color', [1 1 1], 'interp','latex');    
    text(e1(1)-0.1, e1(2)-0.1, e1(3),'$E_1$','fontsize',14,'interp','latex');    
    text(e2(1)-0.1, e2(2)-0.1, e2(3),'$E_2$','fontsize',14,'interp','latex');    
    text(e2s(1)-0.1, e2s(2)-0.1, e2s(3),'$\tau_{1/4}E_2$','fontsize',14,'interp','latex');    
    text(e3(1)-0.1, e3(2)-0.1, e3(3),'$E_3$','fontsize',14,'interp','latex');
% %     if 0, phase = 0:0.01:1;  ve = eq(3).a(1:2:end)+1i*eq(3).a(2:2:end);
% %     ve = repmat(ve,1,length(phase)).*exp(-2i*pi*(1:31)'*phase);
% %     ae = zeros(size(aa,1),length(phase));
% %     ae(1:2:end,:) = real(ve); ae(2:2:end,:) = imag(ve); ve = v'*ae;
% %     plot3(ve(1,:),ve(2,:),ve(3,:),'k-','linewidth',1.5); end    
    xlabel('$v_1$','fontsize',14,'interp','latex');
    ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
    zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);
%     xlim([-0.1 1.6]);
%     ylim([-0.4 1.2]);
%     zlim([-0.02 0.08]);
     view(86.5, 24);
% %     set(gca,'FontSize',14)
%     arr1 = annotation(1,'arrow',[0.458 0.5347],[0.7372 0.5791]);
%     print -depsc2 ks22_TW1_manifold_mf1.eps

  
    
 %% Add rpos to previous plot.

load ks22f90h25.mat; np=1;h = 0.1; 
refpo=1;
a0=rpo(refpo).a1;
[tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
aa0inv=mfTraj(aa0,slc);
av0=vInv'*aa0inv;
plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [218/255 112/255 214/255 ]);
% refpo=3;
% a0=rpo(refpo).a1;
% [tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
% aa0inv=mfTraj(aa0,slc);
% av0=vInv'*aa0inv;
% plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [189/255 183/255 107/255]);
% refpo=4;
% a0=rpo(refpo).a1;
% [tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
% aa0inv=mfTraj(aa0,slc);
% av0=vInv'*aa0inv;
% plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [139/255 69/255 19/255]);
% refpo=5;
% a0=rpo(refpo).a1;
% [tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
% aa0inv=mfTraj(aa0,slc);
% av0=vInv'*aa0inv;
% plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [255/255 165/255 0/255]);
print -depsc2 ks22_TW1_manifold_rpo1_mf1.eps




%% Add plot of E2 unstable manifold to previous figure
  k = 2;  h = 0.1;  tend = 100;  av = [];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v2 = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1))  real(eq(k).evec(:,7))]);
  for delta = [0:0.03:1, 0.74 0.066 0.0693]*ere*period,
    a0 = eq(k).a + 1e-4.*exp(delta).*v2(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
%    [aa, ss] = ksfm_so2(aa, L);
    [aa, ss] = mfTraj(aa,slc);
    nj = 100;  tj = zeros(1,nj);
    ij = find(abs(diff(ss))>0.9*pi/slc);
    tj(1,1:length(ij)) = (tt(ij)+tt(ij+1))/2;
    tas = [tj tt];  tas = [tas; [NaN(size(aa,1),nj) aa]]';
    tas = sortrows(tas,1)';  aa = tas(2:end,:);
    av = [av; vInv'*aa];
  end
plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','-','color', [	135/255 206/255 250/255]);   

print -depsc2 ks22_TW1_manifold_rpo1_E2_manif_mf1.eps


%% Plot homoclinic rpo in c_i=0, b_i>0 slice
clear; load kse22orbits;  k = 2;  h = 0.1;
v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
slc=1;    
vInv=mfTraj([eq(k).a eq(k).a eq(k).a]+v,slc);
figure(3); clf; view(14,18);
hold on; grid on; % view(27,30);
r1 = vInv'*mfTraj(re(1).a,slc); plot3(r1(1),r1(2),r1(3),'k.','markersize',14);
[e1 s] = mfTraj(eq(1).a,slc);
e1 = vInv'*e1; plot3(e1(1),e1(2),e1(3),'k.','markersize',14);
e2 = vInv'*fmRot(eq(2).a,s); plot3(e2(1),e2(2),e2(3),'k.','markersize',14);
e2s = vInv'*fmRot(eq(2).a,s+pi/2); plot3(e2s(1),e2s(2),e2s(3),'k.','markersize',14);
e3 = vInv'*fmRot(eq(3).a,s); plot3(e3(1),e3(2),e3(3),'k.','markersize',14);
text(r1(1)-0.05, r1(2)-0.05, r1(3)-0.05,'$TW_1$','fontsize',14,  'interp','latex');    
text(e1(1)-0.05, e1(2)-0.05, e1(3)-0.05,'$E_1$','fontsize',14,'interp','latex');    
text(e2(1)-0.05, e2(2)-0.05, e2(3)-0.05,'$E_2$','fontsize',14,'interp','latex');    
text(e2s(1)-0.05, e2s(2)-0.05, e2s(3)-0.05,'$\tau_{1/4}E_2$','fontsize',14,'interp','latex');    
text(e3(1)-0.05, e3(2)-0.05, e3(3)-0.05,'$E_3$','fontsize',14,'interp','latex');
xlabel('$v_1$','fontsize',14,'interp','latex');
ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);


load ks22f90h25.mat; np=1;h = 0.01; 

refpo=62;
a0=rpo(refpo).a1;
[tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
aa0inv=mfTraj(aa0,slc);
av0=vInv'*aa0inv;
plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [255/255 165/255 0/255]);

refpo=5;
a0=rpo(refpo).a1;
[tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
aa0inv=mfTraj(aa0,slc);
av0=vInv'*aa0inv;
plot3(av0(1,:)',av0(2,:)',av0(3,:)','r.-');

%% Add plot of E2 unstable manifold to previous figure
  k = 2;  h = 0.01;  tend = 100;  av = [];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v2 = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1))  real(eq(k).evec(:,7))]);
  for delta = [0:0.06:1, 0.74 0.066 0.0693]*ere*period,
    a0 = eq(k).a + 1e-4.*exp(delta).*v2(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 1);
%    [aa, ss] = ksfm_so2(aa, L);
    [aa, ss] = mfTraj(aa,slc);
    nj = 100;  tj = zeros(1,nj);
    ij = find(abs(diff(ss))>0.9*pi/slc);
    tj(1,1:length(ij)) = (tt(ij)+tt(ij+1))/2;
    tas = [tj tt];  tas = [tas; [NaN(size(aa,1),nj) aa]]';
    tas = sortrows(tas,1)';  aa = tas(2:end,:);
    av = [av; vInv'*aa];
  end
plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','-','color', [	135/255 206/255 250/255]);   

print -depsc2 ks22_E2_manif_homo_rpo_mf1.eps

%% TW1 unstable manifold using slice c_i-0, b_i>0, Fourier modes projection (based on ks22figs.m) 
  clear; load kse22orbits;  k = 1;  h = 0.1;  tend = 90;  av = [];
  ere = real(re(k).eig(1));  period = 2*pi/imag(re(k).eig(1));
  v = gsorth([real(re(k).evec(:,1)) imag(re(k).evec(:,1)) real(re(k).evec(:,3))]);
  slc=1;
  % Pick projection
  proj=[3,7,8]; 
  for delta = [0:0.03:1, 0.74 0.066 0.0693]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
%    [aa, ss] = ksfm_so2(aa, L);
    [aa, ss] = mfTraj(aa,slc);
    nj = 100;  tj = zeros(1,nj);
    ij = find(abs(diff(ss))>0.9*pi/slc);
    tj(1,1:length(ij)) = (tt(ij)+tt(ij+1))/2;
    tas = [tj tt];  tas = [tas; [NaN(size(aa,1),nj) aa]]';
    tas = sortrows(tas,1)';  aa = tas(2:end,:);
    av = [av; aa(proj(1),:); aa(proj(2),:); aa(proj(3),:);];
  end
  figure(1); clf; %set(gcf,'pos',[5 400 420 300],'paperpos',[5 8 12 9]); 
% %   ax1 = axes('pos',[0.12 0.10 0.83 0.88]);
    plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','-','color',[.5 .5 .5]);
    hold on; grid on; % view(27,30);
    plot3(av(end-8,:)',av(end-7,:)',av(end-6,:)','.-','color',[0 .8 0]);
     plot3(av(end-5,:)',av(end-4,:)',av(end-3,:)','r.-');
%     plot3(av(end-2,:)',av(end-1,:)',av(end,:)','b.-');
% %     text(0.5, 0.1,-0.5,'A','fontsize',14,'color',[0 .8 0]);
% %     text(1.0, 0.0,-0.5,'B','fontsize',14,'color','r');
% %     text(1.0, 0.0, 0.3,'C','fontsize',14,'color','b');
%    text(-0.4, 0.3,-0.15,'$\mathbf{g}(L/4)E_2$','fontsize',14,'interp','latex');    
%     text(-0.6, 0.3, 0.2,'$\tau_{1/4}E_2$','fontsize',14,'interp','latex');    
%     text(1.1, 0.0,-0.1,'$E_3$','fontsize',14,'interp','latex');
%     r1 = mfTraj(re(1).a,slc); plot3(r1(1),r1(2),r1(3),'k.','markersize',14);
%     e1 = mfTraj(eq(1).a,slc); plot3(e1(1),e1(2),e1(3),'k.','markersize',14);
%     e2 = mfTraj(eq(2).a,slc); plot3(e2(1),e2(2),e2(3),'k.','markersize',14);
%     e3 = mfTraj(eq(3).a,slc); plot3(e3(1),e3(2),e3(3),'k.','markersize',14);
%     text(r1(1)-0.1, r1(2)-0.1, r1(3),'$TW_1$','fontsize',14, 'color', [1 1 1], 'interp','latex');    
%     text(e1(1)-0.1, e1(2)-0.1, e1(3),'$E_1$','fontsize',14,'interp','latex');    
%     text(e2(1)-0.1, e2(2)-0.1, e2(3),'$E_2$','fontsize',14,'interp','latex');    
%     text(e3(1)-0.1, e3(2)-0.1, e3(3),'$E_3$','fontsize',14,'interp','latex');
% % %     if 0, phase = 0:0.01:1;  ve = eq(3).a(1:2:end)+1i*eq(3).a(2:2:end);
% %     ve = repmat(ve,1,length(phase)).*exp(-2i*pi*(1:31)'*phase);
% %     ae = zeros(size(aa,1),length(phase));
% %     ae(1:2:end,:) = real(ve); ae(2:2:end,:) = imag(ve); ve = v'*ae;
% %     plot3(ve(1,:),ve(2,:),ve(3,:),'k-','linewidth',1.5); end    
    xlabel('$v_1$','fontsize',14,'interp','latex');
    ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
    zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);
%     xlim([-0.1 1.6]);
%     ylim([-0.4 1.2]);
%     zlim([-0.02 0.08]);
     view( 1.245000000000000e+02,20);
% %     set(gca,'FontSize',14)
%     arr1 = annotation(1,'arrow',[0.458 0.5347],[0.7372 0.5791]);
    print -depsc2 ks22_TW1_manifold_mf1.eps

%% E2 unstable manifold using slice c_i=0, b_i>0, Fourier modes projection (based on ks22figs.m) 
  clear;  load kse22orbits;  k = 2;  h = 0.1;  tend = 150;  av = [];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
  % Pick slice:
  slc=1;
  % Pick projection
  proj=[1,4,5]; 
  for delta = [0:0.03:1, 0.74 0.066 0.0693]*ere*period,
    a0 = eq(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
%    [aa, ss] = ksfm_so2(aa, L);
    [aa, ss] = mfTraj(aa,slc);
     nj = 100;  tj = zeros(1,nj);
     ij = find(abs(diff(ss))>0.9*pi/slc);
     tj(1,1:length(ij)) = (tt(ij)+tt(ij+1))/2;
     tas = [tj tt];  tas = [tas; [NaN(size(aa,1),nj) aa]]';
     tas = sortrows(tas,1)';  aa = tas(2:end,:);
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

%% E2 unstable manifold in invariant coordinates, "Fourier modes" projection (based on ks22figs.m) 
  clear;  load kse22orbits;  k = 2;  h = 0.1;  tend = 150;  av = [];
  v = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1)) real(eq(k).evec(:,7))]);
   proj=[1,3,5];
% proj=[1,4,5];
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
  figure(1);clf;% set(gcf,'pos',[5 400 420 300],'paperpos',[5 8 12 9]); 
% %   ax1 = axes('pos',[0.12 0.10 0.83 0.88]);
    plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','-','color',[.5 .5 .5]);
    hold on; grid on; % view(27,30);
%     plot3(av(end-8,:)',av(end-7,:)',av(end-6,:)','.-','color',[0 .8 0]);
%     plot3(av(end-5,:)',av(end-4,:)',av(end-3,:)','r.-');
%     plot3(av(end-2,:)',av(end-1,:)',av(end,:)','b.-');
% %     text(0.5, 0.1,-0.5,'A','fontsize',14,'color',[0 .8 0]);
% %     text(1.0, 0.0,-0.5,'B','fontsize',14,'color','r');
% %     text(1.0, 0.0, 0.3,'C','fontsize',14,'color','b');
%    text(-0.4, 0.3,-0.15,'$\mathbf{g}(L/4)E_2$','fontsize',14,'interp','latex');    
%     text(-0.6, 0.3, 0.2,'$\tau_{1/4}E_2$','fontsize',14,'interp','latex');    
%     text(1.1, 0.0,-0.1,'$E_3$','fontsize',14,'interp','latex');
    e2 = mfinvSO2(eq(k).a); plot3(e2(proj(1)),e2(proj(2)),e2(proj(3)),'k.','markersize',30);
    e3 = mfinvSO2(eq(3).a); plot3(e3(proj(1)),e3(proj(2)),e3(proj(3)),'k.','markersize',30);
    text(e2(proj(1)), e2(proj(2)), e2(proj(3))-0.1,'$E_2$','fontsize',14,'interp','latex');    
    text(e3(proj(1)), e3(proj(2))-0.05, e3(proj(3)),'$E_3$','fontsize',14,'interp','latex');
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
     view(-88.5,-10);
%     arr1 = annotation(1,'arrow',[0.458 0.5347],[0.7372 0.5791]);
     print -depsc2 ks22_E2_manifold_proj145.eps

%% TW1 unstable manifold in invariant coordinates, "Fourier modes" projection (based on ks22figs.m) 
  clear;  load kse22orbits;  k = 1;  h = 0.1;  tend = 70;  av = [];
  v = gsorth([real(re(k).evec(:,1)) imag(re(k).evec(:,1)) real(re(k).evec(:,2))]);
  proj=[1,5,9];
% proj=[1,4,5];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(re(k).eig(1));
  for delta = [0:0.03:1, 0.74 0.51 0.54]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
    if delta == 0.74*ere*period, aa1 = aa; end
    if delta == 0.066*ere*period, aa2 = aa; end, 
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av = [av; aa(proj(1),:); aa(proj(2),:); aa(proj(3),:);];
  end
  figure(1);clf;% set(gcf,'pos',[5 400 420 300],'paperpos',[5 8 12 9]); 
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
    e2 = mfinvSO2(eq(2).a); plot3(e2(proj(1)),e2(proj(2)),e2(proj(3)),'k.','markersize',30);
    e3 = mfinvSO2(eq(3).a); plot3(e3(proj(1)),e3(proj(2)),e3(proj(3)),'k.','markersize',30);
    text(e2(proj(1)), e2(proj(2)), e2(proj(3))-0.1,'$E_2$','fontsize',14,'interp','latex');    
    text(e3(proj(1)), e3(proj(2))-0.05, e3(proj(3)),'$E_3$','fontsize',14,'interp','latex');
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
%      view(-88.5,-10);
%     arr1 = annotation(1,'arrow',[0.458 0.5347],[0.7372 0.5791]);
     print -depsc2 ks22_TW1_manifold_proj145.eps
     
    figure(2);clf; 
    difs=zeros(size(aa,2),1);
    difs0=zeros(size(aa,2),1);
    for i=1:size(difs,1),
        difs(i)=norm(mfinvSO2(aa1(:,i))-mfinvSO2(re(k).a));        
        difs0(i)=norm(aa1(:,i)-re(k).a);        
    end
    semilogy(difs/max(difs),'.'); hold on;
    semilogy(difs0/max(difs0),'r.');

%% Add plot of E2 unstable manifold in previous figure
  k = 2;  h = 0.1;  tend = 150;  av = [];
  ere = real(eq(k).eig(1));  period = 2*pi/imag(eq(k).eig(1));
  v2 = gsorth([real(eq(k).evec(:,1)) imag(eq(k).evec(:,1))  real(eq(k).evec(:,7))]);
  for delta = [0:0.03:1, 0.74 0.066 0.0693]*ere*period,
    a0 = eq(k).a + 1e-4.*exp(delta).*v2(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 3);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    av = [av; aa(proj(1),:); aa(proj(2),:); aa(proj(3),:);];
    if delta == 0.74*ere*period, aa1 = aa; end
    if delta == 0.066*ere*period, aa2 = aa; end, 
  end
plot3(av(1:3:end-8,:)',av(2:3:end-7,:)',av(3:3:end-6,:)','m-');


%% Shadowing of RPO(64.51) in physical space


load ks22f90h25.mat;
h=0.1; L=22;
np=1;

refpo=36; % Pick reference cycle

% Integrate reference cycle i.c.
a0=rpo(refpo).a1;
[tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 

% Real space plot

[x, uu] = ksfm2real(aa0, L, 64);
figure(); set(gcf,'pos',[100 500 250 350]); clf;
ax1 = axes('pos',[0.22 0.12 0.70 0.78]); pcolor(x,tt0,uu'); caxis([-3 3]);
shading interp; colormap('jet');
xlabel('x','fontsize',14);  ylabel('t','fontsize',14,'rotat',0);
T=rpo(refpo).T1; s=rpo(refpo).s1;
title(['T_p=' num2str(T) ', d_p=' num2str(s) ], 'fontsize', 14)
% set(get(gca,'ylabel'),'pos',[-13.5 74 1]);
% set(gcf,'paperpos',[8 10 6 10]); 


exportfig(gcf, 'ks22rpoT64p51s2p5phys.eps', 'Color', 'rgb', 'resolution', 300);

% print -depsc2  'ks22rpoT64p51s2p5phys.eps';

% Pick shadowing cycle 

ipo=1;

a0 = rpo(ipo).a1;
[tt, aa] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 

% Real space plot

[x, uu] = ksfm2real(aa, L, 64);
figure(); set(gcf,'pos',[100 500 250 350]); clf;
ax1 = axes('pos',[0.22 0.12 0.70 0.78]); pcolor(x,tt,uu'); caxis([-3 3]);
shading interp; colormap('jet'); hold on;
xlabel('x','fontsize',14);  ylabel('t','fontsize',14,'rotat',0);
T=rpo(ipo).T1; s=rpo(ipo).s1;
title(['T_p=' num2str(T) ', d_p=' num2str(s) ], 'fontsize', 14)
ne = ceil(rpo(refpo).T1/T);
plot(x([1 end])*ones(1,ne),[T;T]*(1:ne),'w-','LineWidth',1.2);
% set(get(gca,'ylabel'),'pos',[-13.5 74 1]);
%set(gcf,'paperpos',[8 10 6 10]); 

exportfig(gcf, 'ks22rpoT16p31s2p9phys.eps', 'Color', 'rgb', 'resolution', 300);

%%%%

%%% plot alternative shadowing ppo
ipo1=17;
a0=ksfmRefl(rpo(ipo1).a1);
[tt1, aa1] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 

% Real space plot

[x, uu] = ksfm2real(aa1, L, 64);
figure(); set(gcf,'pos',[100 500 250 350]); clf;
ax1 = axes('pos',[0.22 0.12 0.70 0.78]); pcolor(x,tt1,uu'); caxis([-3 3]);
shading interp; colormap('jet'); hold on;
xlabel('x','fontsize',14);  ylabel('t','fontsize',14,'rotat',0);
T=rpo(ipo1).T1; s=rpo(ipo1).s1;
title(['T_p=' num2str(T) ', d_p=' num2str(s) ],'fontsize',14);
ne = ceil(rpo(refpo).T1/T);
plot(x([1 end])*ones(1,ne),[T;T]*(1:ne),'w-','LineWidth',1.2);
%set(gcf,'paperpos',[8 10 6 10]); 

exportfig(gcf,'ks22rpoT47p32s0p3phys.eps', 'Color', 'rgb', 'resolution', 300);

figure();
set(gcf,'pos',[100 150 250 150]); clf;
% ax1 = axes('pos',[0.22 0.12 0.70 0.78]); 
plot(x,uu(:,1)','LineWidth',2);  xlim([-11,11]); ylim([-2.2 2.2]);
xlabel('x','fontsize',12);  ylabel('u(x,0)','fontsize',12);

exportfig(gcf,'ks22rpoT47p32s0p3templ.eps', 'Color', 'rgb', 'resolution', 300);


%% Shadowing of RPO(64.51)

clear; load kse22orbits;
k=1;
v = gsorth([real(re(k).evec(:,1)) imag(re(k).evec(:,1)) real(re(k).evec(:,3))]);
vInv=mfinvTraj(repmat(re(k).a,1,3)+v);

load ks22f90h25.mat;
h=0.1; L=22;
np=1;

refpo=36; % Pick reference cycle

% Integrate reference cycle i.c.
a0=rpo(refpo).a1;
[tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 
save('ks22rpo77p85s10p8.mat','aa0');
% Pick shadowing cycle 

ipo=1;

a0 = rpo(ipo).a1;
[tt, aa] = ksfmetd2(a0, L, h, rpo(ipo).T1, np); 

% % Compute distance of points on cycles
% [d, aamf, aa0mf, pos]=minDistanceInvPos(aa,aa0);

% globmin=find(d==min(d));

% points of minimal distance (to plot if needed)
% if size(aa0mf,2)>=size(aamf,2), 
%     pnt0mf = aa0mf(:,pos(globmin));
%     pntmf = aamf(:,globmin);
% else
%     pnt0mf = aa0mf(:,globmin);
%     pntmf = aamf(:,pos(globmin));
% end

% % difference vector at minimal distance
% difvmf=pntmf-pnt0mf;

% 
% % pick projection for plot
% proj=[ 1, 27, 28];


%%%%

%%% plot alternative shadowing ppo
ipo1=17;
a0=ksfmRefl(rpo(ipo1).a1);
[tt, aa1] = ksfmetd2(a0, L, h, rpo(ipo1).T1, np); 

save('ks22rpo55p6s5p25.mat','aa1');

[d1, aamf1, aa0mf, pos1]=minDistanceInvPos(aa1,aa0);

% globmin=find(d1==min(d1));
% 
% % points of minimal distance (to plot if needed)
%  if size(aa0mf,2)>=size(aamf1,2), 
%      pnt0mf = aa0mf(:,pos1(globmin));
%      pntmf = aamf(:,globmin);
%  else
%      pnt0mf = aa0mf(:,globmin);
%      pntmf = aamf1(:,pos(globmin));
%  end

% % difference vector at minimal distance
% difvmf=pntmf-pnt0mf;

%%%%
% %%% plot alternative ppo
% ipo2=7;
% a0=ppo(ipo2).a;
% [tt, aa2] = ksfmetd2(a0, L, h, ppo(ipo2).T, np); 
% 
% [d, aamf2, aa0mf]=minDistanceInv(aa2,aa0);
% 
% plot3(aamf2(proj(1),:),aamf2(proj(2),:),aamf2(proj(3),:),'k-', 'LineWidth', 2);
% %%%%
xlabel('$v_1$','fontsize',14,'interp','latex');
ylabel('$v_2$','fontsize',14,'interp','latex');
zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');
set(gca,'FontSize',14);
axis tight;

l2=['T_p=' num2str(rpo(refpo).T1,4)];
l1=['T_p=' num2str(rpo(ipo).T1,4)];
l3=['T_p=' num2str(rpo(ipo1).T1,4)];
% l4=['T_p=' num2str(ppo(ipo2).T,4)];


% use [az, el]=view; 

%%%%%

% Compute Floquet exponents (from stored multipliers)
fexp_refpo = [log(abs(rpo(refpo).e(1)))/rpo(refpo).T, log(abs(rpo(refpo).e(2)))/rpo(refpo).T];
fexp_ipo = [log(abs(rpo(ipo).e(1)))/rpo(ipo).T, log(abs(rpo(ipo).e(2)))/rpo(ipo).T];
fexp_ipo1 = [log(abs(rpo(ipo1).e(1)))/rpo(ipo1).T, log(abs(rpo(ipo1).e(4)))/rpo(ipo1).T];


%%%%%

% Changing input data here!!!
load kse22orbits;
% vInv=gsorth([mfinvTraj(rpo(1).evec(:,1)) vInv(:,1) vInv(:,2)]);
vInv=gsorth([mfinvTraj(rpo(1).evec(:,1)) mfinvTraj(rpo(1).evec(:,4)) vInv(:,1)]);
% vInv=gsorth([vPoinc vInv(:,1) vInv(:,2)]);

aamf=vInv'*mfinvTraj(aa);
aa0mf=vInv'*mfinvTraj(aa0);
aamf1=vInv'*mfinvTraj(aa1);


fig2= figure();

plot3(aamf(1,:),aamf(2,:),aamf(3,:),'s-');

hold on;

plot3(aa0mf(1,:),aa0mf(2,:),aa0mf(3,:),'r.-');


plot3(aamf1(1,:),aamf1(2,:),aamf1(3,:),'g-', 'LineWidth', 2);


axis tight;
xlabel('$v_1$','fontsize',14,'interp','latex');
ylabel('$v_2$','fontsize',14,'interp','latex');
zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');
set(gca,'FontSize',14);

view(-101.5,4)
axis tight;

legend(l1,l2,l3);

print -depsc2 'ks22rpoT6451shad.eps';

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

saveas(fig4, ['ks22rpoT', '6451', 'angl_dist.png']);
saveas(fig4, ['ks22rpoT', '6451', 'angl_dist.pdf']);
saveas(fig4, ['ks22rpoT', '6451', 'angl_dist.eps']);

%% Shadowing of RPO(80.71)

clear; load kse22orbits;
k=1;
ere = real(re(k).eig(1));  period = 2*pi/imag(re(k).eig(1));
v = gsorth([real(re(k).evec(:,1)) imag(re(k).evec(:,1)) real(re(k).evec(:,3))]);
vInv=mfinvTraj(repmat(re(k).a,1,3)+v);

load ks22f90h25.mat;
h=0.1; L=22;
np=1;

refpo=105; % Pick reference cycle

% Integrate reference cycle i.c.
a0=rpo(refpo).a1;
[tt0, aa0] = ksfmetd2(a0, L, h, rpo(refpo).T1, np); 


% Pick shadowing cycle
ipo=1;

a0 = rpo(ipo).a1;
[tt, aa] = ksfmetd2(a0, L, h, rpo(ipo).T1, np); 

% % Compute distance of points on cycles
% [d, aamf, aa0mf, pos]=minDistanceInvPos(aa,aa0);

% globmin=find(d==min(d));

% points of minimal distance (to plot if needed)
% if size(aa0mf,2)>=size(aamf,2), 
%     pnt0mf = aa0mf(:,pos(globmin));
%     pntmf = aamf(:,globmin);
% else
%     pnt0mf = aa0mf(:,globmin);
%     pntmf = aamf(:,pos(globmin));
% end

% % difference vector at minimal distance
% difvmf=pntmf-pnt0mf;

% 
% % pick projection for plot
% proj=[ 1, 27, 28];


%%% alternative shadowing ppo
ipo1=17;
a0=ksfmRefl(rpo(ipo1).a1);
[tt, aa1] = ksfmetd2(a0, L, h, rpo(ipo1).T1, np); 


fexp_refpo = [log(abs(ppo(refpo).e(1)))/ppo(refpo).T, log(abs(ppo(refpo).e(2)))/ppo(refpo).T];
fexp_ipo = [log(abs(rpo(ipo).e(1)))/rpo(ipo).T, log(abs(rpo(ipo).e(2)))/rpo(ipo).T];
fexp_ipo1 = [log(abs(rpo(ipo1).e(1)))/rpo(ipo1).T, log(abs(rpo(ipo1).e(4)))/rpo(ipo1).T];

pPoinc=mfinvTraj(aa(:,2));
vPoinc= mfinvTraj(aa(:,3))-mfinvTraj(aa(:,2));
dir=1;
pPoinc1=mfinvTraj(aa1(:,5));
vPoinc1= mfinvTraj(aa1(:,6))-mfinvTraj(aa1(:,5));


l2=['T_p=' num2str(rpo(refpo).T1,4)];
l1=['T_p=' num2str(rpo(ipo).T1,4)];
l3=['T_p=' num2str(rpo(ipo1).T1,4)];

% Changing input data here!!!
load kse22orbits;
vInv=gsorth([mfinvTraj(rpo(1).evec(:,1)) vInv(:,1) vInv(:,2)]);
vInv=gsorth([mfinvTraj(rpo(1).evec(:,1)) mfinvTraj(rpo(1).evec(:,4)) vInv(:,1)]);
vInv=gsorth([vPoinc vInv(:,1) vInv(:,2)]);

aamf=vInv'*mfinvTraj(aa);
aa0mf=vInv'*mfinvTraj(aa0);
aamf1=vInv'*mfinvTraj(aa1);


fig2= figure();

plot3(aamf(1,:),aamf(2,:),aamf(3,:),'s-');

hold on;

plot3(aa0mf(1,:),aa0mf(2,:),aa0mf(3,:),'r.-');


plot3(aamf1(1,:),aamf1(2,:),aamf1(3,:),'g-', 'LineWidth', 2);

% pPoincv=vInv'*pPoinc;
% plot3(pPoincv(1),pPoincv(2),pPoincv(3),'bs','MarkerSize',12)

%%%%
% %%% plot alternative ppo
% ipo2=7;
% a0=ppo(ipo2).a;
% [tt, aa2] = ksfmetd2(a0, L, h, ppo(ipo2).T, np); 
% 
% [d, aamf2, aa0mf]=minDistanceInv(aa2,aa0);
% 
% plot3(aamf2(proj(1),:),aamf2(proj(2),:),aamf2(proj(3),:),'k-', 'LineWidth', 2);
% %%%%
xlabel('$v_1$','fontsize',14,'interp','latex');
ylabel('$v_2$','fontsize',14,'interp','latex');
zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');
set(gca,'FontSize',14);

view(-101.5,4)
axis tight;


% l4=['T_p=' num2str(ppo(ipo2).T,4)];
legend(l1,l2,l3);


print -depsc2  'ks22rpoT8072shad.eps';
%%%%%

% Compute Floquet exponents (from stored multipliers)


%% Add TW1 unstable manifold
  v = gsorth([real(re(k).evec(:,1)) imag(re(k).evec(:,1)) real(re(k).evec(:,3))]);
  vInv=mfinvTraj(repmat(re(k).a,1,3)+v);



% shadowing shortest rpo
  tend=111; av81=[]; aa81=[];
  for delta = [0.50915:0.00001:0.50923]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    aa81 = [aa81; aa];
    av81 = [av81; vInv'*aa];
  end
  
% shadowing shortest rpo
  tend=98; av82=[]; aa82=[];
  for delta = [0.509:0.0005:0.51]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    aa82 = [aa82; aa];
    av82 = [av82; vInv'*aa];
  end  
  tend=95; av8=[]; aa8=[];
  for delta = [0.495:0.0005:0.5105]*ere*period,
    a0 = re(k).a + 1e-4.*exp(delta).*v(:,2);
    [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%    [aa, ss] = ksfm_so2(aa, L);
    aa = mfinvTraj(aa);
    aa8 = [aa8; aa];
    av8 = [av8; vInv'*aa];
  end
% 

save('ks22tw1manifS1.mat','aa81');
save('ks22tw1manifS2.mat','aa82');
save('ks22tw1manifS3.mat','aa8');


plot3(av82(1:3:end,:)',av82(2:3:end,:)',av82(3:3:end,:)','k-'); % shadow rpo 1
 plot3(av81(1:3:end,:)',av81(2:3:end,:)',av81(3:3:end,:)','k-'); %, 'color', [0.5 0.5 0.5]); % shadow rpo 1
xlabel('$v_1$','fontsize',14,'interp','latex');
ylabel('$v_2$','fontsize',14,'interp','latex');%,'pos',[1.3 0 -.8]); 
zlabel('$v_3$','rotat',0,'fontsize',14,'interp','latex');%,'pos',[-1.4 -.7 0]);
set(gca,'FontSize',14)
% load ks22f90h25.mat; np=1; h = 0.1; 
% refpo=1;
% a0=rpo(refpo).a1;
% [tt0, aa0] = ksfmetd(a0, L, h, rpo(refpo).T1, np); 
% aa0inv=mfinvTraj(aa0);
% av0=vInv'*aa0inv;
% plot3(av0(1,:)',av0(2,:)',av0(3,:)','.-', 'color', [218/255 112/255 214/255 ]);
print -depsc2 ks22_TW1_manif_rpos_inv.eps
% plot3(aamf(1,1),aamf(2,1),aamf(3,1),'k*');
% plot3(aamf(1,2),aamf(2,2),aamf(3,2),'kv');

% Add RPO(16.31) unstable manifold
% load kse22orbits;
ere = fexp_ipo(1);  period = rpo(1).T;
tend=20.*period; av81=[]; av81p=[]; av81p2=[]; a81p=[]; a81p2=[];
nhit=2;
for delta = [0.:0.005:0.999]*ere*period,
    a0 = rpo(1).a - exp(delta)*1e-4.*rpo(1).evec(:,1);
    [aPoinc, aa] = ksfmIntPoinc(a0, pPoinc, vPoinc, dir, L, h, nhit, tend);
    a81p=[a81p, vInv'*aPoinc];
    av81=[av81, vInv'*aa, [NaN, NaN, NaN]'];
    av81p = [av81p, reshape(vInv'*aPoinc,[],1)];
%     [aPoinc, aa] = ksfmIntPoinc(a0, pPoinc1, vPoinc1, dir, L, h, nhit, tend);
%     a81p2=[a81p2, vInv'*aPoinc];
%     av81p2 = [av81p2, reshape(vInv'*aPoinc,[],1)];
end
a0 = rpo(1).a;
[aPoinc, aa] = ksfmIntPoinc(a0, pPoinc, vPoinc, dir, L, h, nhit, tend);
ac1p= vInv'*aPoinc;
avc1=vInv'*aa;
avc1p = reshape(vInv'*aPoinc,[],1);
[a81ps1, a81ps1ds]=keepCont(a81p(:,1:nhit:end),0.1);
[a81ps2, a81ps2ds]=keepCont(a81p(:,2:nhit:end),0.1);
% plot3(av81(1:3:end,:)',av81(2:3:end,:)',av81(3:3:end,:)','c-'); 
% plot3(av81p(1:3:end,:)',av81p(2:3:end,:)',av81p(3:3:end,:)','c*','MarkerSize',10); 
% plot3(av81p2(1:3:end,:)',av81p2(2:3:end,:)',av81p2(3:3:end,:)','cs','MarkerSize',10); 
plot3(a81ps1(1,:)',a81ps1(2,:)',a81ps1(3,:)','c.','MarkerSize',10); 
plot3(a81ps2(1,:)',a81ps2(2,:)',a81ps2(3,:)','m.','MarkerSize',10); 
% plot3(a81p(1,4:nhit:size(a81ps2,2))',a81p(2,4:nhit:size(a81ps2,2))',a81p(3,4:nhit:size(a81ps2,2))','rx','MarkerSize',12); 
% plot3(a81p(1,4:nhit:end)',a81p(2,4:nhit:end)',a81p(3,4:nhit:end)','bx','MarkerSize',12); 
% plot3(a81p2(1,1:nhit:end)',a81p2(2,1:nhit:end)',a81p2(3,1:nhit:end)','c.','MarkerSize',10); 
% plot3(a81p(1,2:nhit:end)',a81p(2,2:nhit:end)',a81p(3,2:nhit:end)','mo','MarkerSize',10); 
% plot3(a81p(1,3:nhit:end)',a81p(2,3:nhit:end)',a81p(3,3:nhit:end)','b*','MarkerSize',12); 
% plot3(a81p(1,4:nhit:end)',a81p(2,4:nhit:end)',a81p(3,4:nhit:end)','r*','MarkerSize',12); 
[rm81p, s81p]=manif2rm(ac1p(:,1),a81p,nhit,1);


tend=10.*period; av82=[]; av82p=[]; av82p2=[]; a82p=[]; a82p2=[];
nhit=2;
for delta = [0.:0.005:0.999]*ere*period,
    a0 = rpo(1).a + exp(delta)*1e-4.*rpo(1).evec(:,1);
    [aPoinc, aa] = ksfmIntPoinc(a0, pPoinc, vPoinc, dir, L, h, nhit, tend);
    a82p=[a82p, vInv'*aPoinc];
    av82=[av82, vInv'*aa, [NaN, NaN, NaN]'];
    av82p = [av82p, reshape(vInv'*aPoinc,[],1)];
%     [aPoinc, aa] = ksfmIntPoinc(a0, pPoinc1, vPoinc1, dir, L, h, nhit, tend);
%     a82p2=[a82p2, vInv'*aPoinc];
%     av82p2 = [av82p2, reshape(vInv'*aPoinc,[],1)];
end
[a82ps1, a82ps1ds]=keepCont(a82p(:,1:nhit:end),0.1);
[a82ps2, a82ps2ds]=keepCont(a82p(:,2:nhit:end),0.1);
% plot3(av82(1:3:end,:)',av82(2:3:end,:)',av82(3:3:end,:)','c-'); 
% plot3(av82p(1:3:end,:)',av82p(2:3:end,:)',av82p(3:3:end,:)','c.','MarkerSize',14); 
% plot3(av82p2(1:3:end,:)',av82p2(2:3:end,:)',av82p2(3:3:end,:)','cs','MarkerSize',10); 
plot3(a82ps1(1,:)',a82ps1(2,:)',a82ps1(3,:)','c.','MarkerSize',10); 
plot3(a82ps2(1,:)',a82ps2(2,:)',a82ps2(3,:)','m.','MarkerSize',10); 
% plot3(a82p(1,3:nhit:end)',a82p(2,3:nhit:end)',a82p(3,3:nhit:end)','b*','MarkerSize',12); 
% plot3(a82p2(1,1:nhit:end)',a82p2(2,1:nhit:end)',a82p2(3,1:nhit:end)','c.','MarkerSize',10); 
% plot3(a82p2(1,2:nhit:end)',a82p2(2,2:nhit:end)',a82p2(3,2:nhit:end)','m.','MarkerSize',10); 
% plot3(a82p2(1,3:nhit:end)',a82p2(2,3:nhit:end)',a82p2(3,3:nhit:end)','b*','MarkerSize',12); 

% tend=3*period; av83=[]; 
% for delta = [1.1:0.1:2.5],
%     a0 = rpo(1).a - delta*1e-4.*rpo(1).evec(:,1)*ere*period;
%     [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%     aa = mfinvTraj(aa);
%     av83 = [av83; vInv'*aa];
% end
% plot3(av83(1:3:end,:)',av83(2:3:end,:)',av83(3:3:end,:)','y-'); 
% 
% tend=4*period; av83=[]; 
% for delta = [2.11:0.01:2.15],
%     a0 = rpo(1).a - delta*1e-4.*rpo(1).evec(:,1)*ere*period;
%     [tt, aa] = ksfmetd2(a0, L, h, tend, 2);
%     aa = mfinvTraj(aa);
%     av83 = [av83; vInv'*aa];
% end
% plot3(av83(1:3:end,:)',av83(2:3:end,:)',av83(3:3:end,:)','y.-'); 

view(-101.5,4)

axis tight;

% Add RPO(47.32) unstable manifold
ere = fexp_ipo1(1);  period = rpo(14).T;
tend=10*period; av91=[]; a91p=[]; av91p=[]; a91p2=[]; av91p2=[]; nhit=4;
for delta = [0.:0.005:0.999]*ere*period,
    a0 = ksfmRefl(rpo(14).a - exp(delta)*1e-4.*rpo(14).evec(:,1));
    [aPoinc, aa] = ksfmIntPoinc(a0, pPoinc, vPoinc, dir, L, h, nhit, tend);
    a91p=[a91p, vInv'*aPoinc];
    av91=[av91, vInv'*aa, [NaN, NaN, NaN]'];
    av91p = [av91p, reshape(vInv'*aPoinc,[],1)];
%     [aPoinc, aa] = ksfmIntPoinc(a0, pPoinc1, vPoinc1, dir, L, h, nhit, tend);
%     a91p2=[a91p2, vInv'*aPoinc];
%     av91p2 = [av91p2, reshape(vInv'*aPoinc,[],1)];
end
[a91ps1, a91ps1ds]=keepCont(a91p(:,1:nhit:end),0.1);
[a91ps2, a91ps2ds]=keepCont(a91p(:,2:nhit:end),0.1);
[a91ps3, a91ps3ds]=keepCont(a91p(:,3:nhit:end),0.1);
[a91ps4, a91ps4ds]=keepCont(a91p(:,4:nhit:end),0.1);
% plot3(av91(1:3:end,:)',av91(2:3:end,:)',av91(3:3:end,:)','m.-'); 
% plot3(a91ps1(1,:)',a91ps1(2,:)',a91ps1(3,:)','b.','MarkerSize',10); 
plot3(a91ps2(1,:)',a91ps2(2,:)',a91ps2(3,:)','b.','MarkerSize',10); 
plot3(a91ps4(1,:)',a91ps4(2,:)',a91ps4(3,:)','k.','MarkerSize',10); 
% plot3(a91ps3(1,:)',a91ps3(2,:)',a91ps3(3,:)','b.','MarkerSize',10); 


tend=10*period; av92=[]; av92p=[]; av92p2=[]; a92p=[]; a92p2=[]; nhit=4;
for delta = [0.:0.005:0.999]*ere*period,
    a0 = ksfmRefl(rpo(14).a + exp(delta)*1e-4.*rpo(14).evec(:,1));
    [aPoinc, aa] = ksfmIntPoinc(a0, pPoinc, vPoinc, dir, L, h, nhit, tend);
    a92p=[a92p, vInv'*aPoinc];
    av92=[av92, vInv'*aa, [NaN, NaN, NaN]'];
    av92p = [av92p, reshape(vInv'*aPoinc,[],1)];
%     [aPoinc, aa] = ksfmIntPoinc(a0, pPoinc1, vPoinc1, dir, L, h, nhit, tend);
%     a92p2=[a92p2, vInv'*aPoinc];
%     av92p2 = [av92p2, reshape(vInv'*aPoinc,[],1)];
end
[a92ps1, a92ps1ds]=keepCont(a92p(:,1:nhit:end),0.1);
[a92ps2, a92ps2ds]=keepCont(a92p(:,2:nhit:end),0.1);
[a92ps3, a92ps3ds]=keepCont(a92p(:,3:nhit:end),0.1);
[a92ps4, a92ps4ds]=keepCont(a92p(:,4:nhit:end),0.1);
% plot3(av92(1:3:end,:)',av92(2:3:end,:)',av92(3:3:end,:)','m-'); 
plot3(a92ps2(1,:)',a92ps2(2,:)',a92ps2(3,:)','b.','MarkerSize',10); 
plot3(a92ps4(1,:)',a92ps4(2,:)',a92ps4(3,:)','k.','MarkerSize',10); 
% plot3(a92ps3(1,:)',a92ps3(2,:)',a92ps3(3,:)','k.','MarkerSize',10); 

print -depsc2  'ks22rpoT8072shad_manif.eps';

% Plot return map
figure(); hold on;
plot(rm81p(:,1),rm81p(:,2),'g.-')
rm81p2=manifRep2rm(a81p(:,4:nhit:end), ac1p(:,1), [a81ps1 a81ps2], s81p, size(a81ps1,2), 0.02);


