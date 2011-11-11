% Visualization of KS L=22 reduced state space.
% ES 2011-11-11

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
[tt, aa0] = ksfmetd(a0, L, h, ppo(refpo).T, np); 

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

[d1, aamf1, aa0mf]=minDistanceInv(aa1,aa0);

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

view(-1.5,-2.)

saveas(fig2, ['ks22rpoT', '7035', 'shad.png']);
%%%%%

% Compute Floquet exponents (from stored multipliers)
fexp_refpo = [log(abs(ppo(refpo).e(1)))/rpo(refpo).T, log(abs(ppo(refpo).e(2)))/ppo(refpo).T];
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

fig3= figure();

plot3(aa0mf(proj(1),:),aa0mf(proj(2),:),aa0mf(proj(3),:),'.-');

hold on;

dmin=min(d);
dmax=max(d);

for i=1:size(aamf,2),
    plot3(aamf(proj(1),i),aamf(proj(2),i),aamf(proj(3),i), '.', 'Color', [d(i)/dmax 0 0]);
end


d1min=min(d1);
d1max=max(d1);

for i=1:size(aamf1,2),
    plot3(aamf1(proj(1),i),aamf1(proj(2),i),aamf1(proj(3),i), '.', 'Color', [0 d1(i)/d1max 0]);
end



plot3(aa0mf(proj(1),minanglpos),aa0mf(proj(2),minanglpos),aa0mf(proj(3),minanglpos),'ks');


