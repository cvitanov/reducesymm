clear;

load ks22f90h25angl.mat

h=0.25; N=16; L=22;

np=1;

refpo=48;

sfile=['data/ks22ppo' num2str(refpo) '_rpo_cr.dat'];

crdat=load(sfile);

fig1= figure();

semilogy(crdat(:,1),crdat(:,2),'.');

a0=ppo(refpo).a;
[tt, aa0] = ksfmetd(a0, L, h, ppo(refpo).T, np); 
% [tt, aa0r] = ksfmetd(ksfmRefl(a0), L, h, ppo(refpo).T, np); 
% aa0=[aa0, aa0r];

% ppos close to refpo=4 : 37; 274; 369; 507; 2279;
% maybe not close to refpo=4 : 37;
% not close to refpo=4 : 144;

% ppos close to refpo=1 : 420;
% ppos maybe close to refpo=1 : ipo=60; 166; 170; 201; 282; 288; 17656; 2273; 10365;
%

% rpos maybe close to refpo=1 : ipo=10568;

% ppos close to refpo=5: ipo=1567;

% ppos close to refpo=7: ipo=48;

% rpos close to refpo=115: ipo=410;

% rpos close to refpo=62: ipo=392;
% ppos close to refpo=62: ipo=18;


ipo=6;

a0 = rpo(ipo).a;
%%% fix the phase to the one of first orbit
% phi = angle(ppo(refpo).a(1)+1i*ppo(refpo).a(2));
% ac=a0(1:2:end)+1i*a0(2:2:end);
% as = exp(-1i*(1:size(ac))'*phi).*ac;
% a0(1:2:end)=real(as);
% a0(2:2:end)=imag(as);
[tt, aa] = ksfmetd(a0, L, h, rpo(ipo).T, np); 

[d, aamf, aa0mf, pos]=minDistanceInvPos(aa,aa0);

globmin=find(d==min(d));

% points of minimal distance
if size(aa0mf,2)>=size(aamf,2), 
    pnt0mf = aa0mf(:,pos(globmin));
    pntmf = aamf(:,globmin);
else
    pnt0mf = aa0mf(:,globmin);
    pntmf = aamf(:,pos(globmin));
end

difvmf=pntmf-pnt0mf;

proj=[ 1, 3 ,4];

fig2= figure();

plot3(aamf(proj(1),:),aamf(proj(2),:),aamf(proj(3),:),'.-');

hold on;

plot3(aa0mf(proj(1),:),aa0mf(proj(2),:),aa0mf(proj(3),:),'r.-');
% plot3(pnt0mf( proj(1)), pnt0mf(proj(2)), pnt0mf(proj(3)),'ks');
% plot3(pntmf( proj(1)), pntmf(proj(2)), pntmf(proj(3)),'ko');
%%%%

%%% test alternative ppo
ipo1=9;
a0=rpo(ipo1).a;
[tt, aa1] = ksfmetd(a0, L, h, rpo(ipo1).T, np); 

[d1, aamf1, aa0mf]=minDistanceInv(aa1,aa0);

plot3(aamf1(proj(1),:),aamf1(proj(2),:),aamf1(proj(3),:),'g-', 'LineWidth', 2);
%%%%
% %%% test alternative ppo
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

l2=['T_p=' num2str(rpo(refpo).T,4)];
l1=['T_p=' num2str(rpo(ipo).T,4)];
l3=['T_p=' num2str(rpo(ipo1).T,4)];
% l4=['T_p=' num2str(ppo(ipo2).T,4)];
legend(l1,l2,l3);

% use [az, el]=view; 

view(-1.5,-2.)

saveas(fig2, 'ks22ppo_rpo_shad1.png');
%%%%%

fig3= figure();

plot3(aa(proj(1),:),aa(proj(2),:),aa(proj(3),:));

hold on;

plot3(aa0(proj(1),:),aa0(proj(2),:),aa0(proj(3),:),'r');

fig4= figure();

[d1, aamf1, aa0mf1, pos1]=minDistanceMFpos(aa,aa0); 

plot3(aamf1(proj(1),:),aamf1(proj(2),:),aamf1(proj(3),:),'.-');

globmin1=find(d1==min(d1));
% points of minimal distance
if size(aa0mf1,2)>=size(aamf1,2), 
    pnt0mf1 = aa0mf1(:,pos1(globmin1));
    pntmf1 = aamf1(:,globmin1);
    pnt0mf1o = aa0mf1(:,pos(globmin));
    pntmf1o = aamf1(:,globmin);
else
    pnt0mf1 = aa0mf1(:,globmin1);
    pntmf1 = aamf1(:,pos1(globmin1));
    pnt0mf1o = aa0mf1(:,globmin);
    pntmf1o = aamf1(:,pos(globmin));
end

difvmf1=(pntmf1-pnt0mf1)/norm(pntmf1-pnt0mf1);
difvmf1o=(pntmf1o-pnt0mf1o)/norm(pntmf1o-pnt0mf1o);

hold on;

plot3(aa0mf1(proj(1),:),aa0mf1(proj(2),:),aa0mf1(proj(3),:),'r.-');

plot3(pnt0mf1( proj(1)), pnt0mf1(proj(2)), pnt0mf1(proj(3)),'ks');
plot3(pntmf1( proj(1)), pntmf1(proj(2)), pntmf1(proj(3)),'ko');
plot3(pnt0mf1o( proj(1)), pnt0mf1o(proj(2)), pnt0mf1o(proj(3)),'gs');
plot3(pntmf1o( proj(1)), pntmf1o(proj(2)), pntmf1o(proj(3)),'go');

%%% Compute Jacobian at point of minimal approach
a0=mf(ppo(refpo).a,1);
a0r=mf(ksfmRefl(a0),1);

[tt, aa, da] = ksfmetd(a0, L, h, ppo(refpo).T, np);
[tt2, aa2, da2] = ksfmetd(a0r, L, h, ppo(refpo).T, np);
[ev0, eils] = eig(da2*da);
eils=diag(eils);

%%%



fig5= figure();

proj=[3, 5, 6];


plot3(aamf1(proj(1),:),aamf1(proj(2),:),aamf1(proj(3),:),'.-');

hold on;

plot3(aa0mf1(proj(1),:),aa0mf1(proj(2),:),aa0mf1(proj(3),:),'r.-');

plot3(pnt0mf1( proj(1)), pnt0mf1(proj(2)), pnt0mf1(proj(3)),'ks');
plot3(pntmf1( proj(1)), pntmf1(proj(2)), pntmf1(proj(3)),'ko');
plot3(pnt0mf1o( proj(1)), pnt0mf1o(proj(2)), pnt0mf1o(proj(3)),'gs');
plot3(pntmf1o( proj(1)), pntmf1o(proj(2)), pntmf1o(proj(3)),'go');


