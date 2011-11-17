clear;

load ks22f90h25.mat


h=0.25; N=16; L=22;

np=1;

refpo=97;

sfile=['data/ks22ppo' num2str(refpo) '_rpo_cr.dat'];

a0=ppo(refpo).a;
[tt, aa0] = ksfmetd(a0, L, h, ppo(refpo).T, np);
dtab=ones(500,4);

for ipo=1:500,
    disp(['Comparing ppo ' num2str(refpo) ' , with rpo ' num2str(ipo)]);
    a0=rpo(ipo).a;
    [tt, aa] = ksfmetd(a0, L, h, rpo(ipo).T, np);
    [d, aamf, aa0mf]=minDistanceInv(aa,aa0);
    dtab(ipo,:)=[ipo, mean(d), min(d), max(d)];
end


save(sfile, 'dtab', '-ascii','-double','-tabs', '-append');



% i=59; j=34; % choose snapshots on 2 periodic orbits
% k=5;q=6; % projection
% 
% Ndots=100;
% crcl=(0:Ndots-1)*2*pi/Ndots;
% 
% figure();
% plot(aa0(k,j)*cos(crcl)-aa0(q,j)*sin(crcl),aa0(k,j)*sin(crcl)+aa0(q,j)*cos(crcl));
% hold on;
% plot(aa(k,i)*cos(crcl)-aa(q,i)*sin(crcl),aa(k,i)*sin(crcl)+aa(q,i)*cos(crcl),'r');
% 
% plot(aa0mf(k,j),aa0mf(q,j),'x');
% plot(aamf(k,i),aamf(q,i),'ro');
% 
% figure();
% 
% plot(aa0mf(k,:),aa0mf(q,:));
% hold on;
% plot(aamf(k,:),aamf(q,:),'r');
