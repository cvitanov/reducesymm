clear;

load ks22f90h25angl.mat

h=0.25; N=16; L=22;

np=1;

refpo=37;

sfile=['data/ks22ppo' num2str(refpo) 'cr.dat'];

crdat=load(sfile);

fig1= figure();

semilogy(crdat(:,1),crdat(:,2),'.');

a0=ppo(refpo).a;
[tt, aa0] = ksfmetd(a0, L, h, ppo(refpo).T, np); 

% ppos close to refpo=4 : 37; 274; 369; 507; 2279;
% maybe not close to refpo=4 : 37;
% not close to refpo=4 : 144;

% ppos close to refpo=1 : 420;
% ppos maybe close to refpo=1 : ipo=60; 166; 170; 201; 282; 288; 17656; 2273; 10365;
%

% rpos maybe close to refpo=1 : ipo=10568;

% ppos close to refpo=5: ipo=1567;
ipo=29;

a0=ppo(ipo).a;
[tt, aa] = ksfmetd(a0, L, h, ppo(ipo).T, np); 

[d, aamf, aa0mf]=minDistanceInv(aa,aa0);

proj=[1, 3, 4];

fig2= figure();

plot3(aamf(proj(1),:),aamf(proj(2),:),aamf(proj(3),:),'.-');

hold on;

plot3(aa0mf(proj(1),:),aa0mf(proj(2),:),aa0mf(proj(3),:),'r.-');
% plot3(aa0mf(proj(1),13),aa0mf(proj(2),13),aa0mf(proj(3),13),'ks');
%%%%

xlabel('\beta_1');
ylabel('\beta_2');
zlabel('\gamma_2');

l1=['T_p=' num2str(ppo(refpo).T)];
l2=['T_p=' num2str(ppo(ipo).T)];
% l3=['T_p=' num2str(rpo(ipo1).T)];
legend(l1,l2);

% use [az, el]=view; 

view(-1.5,-2.)

saveas(fig2, 'ks22ppo1_shad.png');
%%%%%

fig3= figure();

plot3(aa(proj(1),:),aa(proj(2),:),aa(proj(3),:));

hold on;

plot3(aa0(proj(1),:),aa0(proj(2),:),aa0(proj(3),:),'r');

[d1, aamf1, aa0mf1]=minDistance(aa,aa0,1); 

fig4= figure();

plot3(aamf1(proj(1),:),aamf1(proj(2),:),aamf1(proj(3),:));

hold on;

plot3(aa0mf1(proj(1),:),aa0mf1(proj(2),:),aa0mf1(proj(3),:),'r');

fig5= figure();

proj=[3, 5, 6];


plot3(aamf1(proj(1),:),aamf1(proj(2),:),aamf1(proj(3),:),'.-');

hold on;

plot3(aa0mf1(proj(1),:),aa0mf1(proj(2),:),aa0mf1(proj(3),:),'r.-');



