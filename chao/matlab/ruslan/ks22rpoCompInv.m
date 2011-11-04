clear;

load ks22f90h25.mat


h=0.25; N=16; L=22;

np=1;

refpo=36;

sfile=['data/ks22rpo' num2str(refpo) 'cr.dat'];

crdat=load(sfile);

fig1= figure();

semilogy(crdat(:,1),crdat(:,2),'.');

a0=rpo(refpo).a;
[tt, aa0] = ksfmetd(a0, L, h, rpo(refpo).T, np); 

% rpos which come close to refpo = 1 : 1330; 105;36; 77
% rpo 36 is shadowed by rpo 1 and rpo 17. ( s = 01)
% rpo 35 is (?) shadowed by rpo 1 and rpo 17. ( s = 01) How to tell apart
% from rpo 36 ?
% rpo 105 is shadowed by rpo 1 (twice) and rpo 17. ( s = 001)
% rpo 6 is close to rpo 47; 309;
% rpo 309 is shadowed by rpo 4 and rpo 6 (twice) (s = '4'6'6')
ipo=1;

a0=rpo(ipo).a;
[tt, aa] = ksfmetd(a0, L, h, rpo(ipo).T, np); 

[d, aamf, aa0mf]=minDistanceInv(aa,aa0);

proj=[ 1, 3, 4];

fig2= figure();

plot3(aa0mf(proj(1),:),aa0mf(proj(2),:),aa0mf(proj(3),:),'r.-');

hold on;

plot3(aamf(proj(1),:),aamf(proj(2),:),aamf(proj(3),:),'bs-');

%%% test alternative rpo
ipo1=17;
a0=rpo(ipo1).a;
[tt, aa1] = ksfmetd(a0, L, h, rpo(ipo1).T, np); 

[d, aamf1, aa0mf]=minDistanceInv(aa1,aa0);

plot3(aamf1(proj(1),:),aamf1(proj(2),:),aamf1(proj(3),:),'g-', 'LineWidth', 2);
%%%%

xlabel('\beta_1');
ylabel('\beta_2');
zlabel('\gamma_2');

l1=['T_p=' num2str(rpo(refpo).T)];
l2=['T_p=' num2str(rpo(ipo).T)];
l3=['T_p=' num2str(rpo(ipo1).T)];
legend(l1,l2,l3);

% use [az, el]=view; 

view(-1.5,-2.)

saveas(fig2, 'ks22rpo_shad1.png');
%%%%%

proj=[ 1, 3, 4];

fig3= figure();

plot3(aa(proj(1),:),aa(proj(2),:),aa(proj(3),:));

hold on;

plot3(aa0(proj(1),:),aa0(proj(2),:),aa0(proj(3),:),'r');

[d1, aamfs, aa0mfs]=minDistanceMF(aa,aa0); 

%%% test alternative rpo
ipo1=17;
a0=rpo(ipo1).a;
[tt, aa1] = ksfmetd(a0, L, h, rpo(ipo1).T, np); 

[d, aamf1s, aa0mfs]=minDistanceMF(aa1,aa0);

%%%%


fig4= figure();

plot3(aamfs(proj(1),:),aamfs(proj(2),:),aamfs(proj(3),:),'.-');

hold on;

plot3(aamf1s(proj(1),:),aamf1s(proj(2),:),aamf1s(proj(3),:),'g.-', 'LineWidth', 2);

plot3(aa0mfs(proj(1),:),aa0mfs(proj(2),:),aa0mfs(proj(3),:),'r.-');

fig5= figure();

proj=[3, 5, 6];

plot3(aamfs(proj(1),:),aamfs(proj(2),:),aamfs(proj(3),:),'.-');

hold on;

plot3(aamf1s(proj(1),:),aamf1s(proj(2),:),aamf1s(proj(3),:),'g.-', 'LineWidth', 2);

plot3(aa0mfs(proj(1),:),aa0mfs(proj(2),:),aa0mfs(proj(3),:),'r.-');




