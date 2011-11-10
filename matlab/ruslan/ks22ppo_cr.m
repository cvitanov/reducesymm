clear;

% load ks22f90h25angl.mat

ph=0; % periodic orbits, only

h=0.25; N=16; L=22;

np=1;

global PFLG;  PFLG = -1; % not needed here?

a0=ppo(1).a;
[tt, aa0] = ksfmetd(a0, L, h, ppo(1).T, np);

for slice=3:3,

for ipo=2:2%size(ppo,2),
    if not(ppo(ipo).angl_prob),
        if ipo~=3 && ipo ~=19. % exclude known repeats (delete from dataset?)
                    disp(['iteration ' num2str(ipo)]);
                    a0=ppo(ipo).a;
                    [tt, aa] = ksfmetd(a0, L, h, ppo(ipo).T, np);
                    [d, aamf, aa0mf]=minDistance(aa,aa0,slice);
                    disp(d);
        end
    end
end

i=59; j=34; % choose snapshots on 2 periodic orbits
k=5;q=6; % projection

Ndots=100;
crcl=(0:Ndots-1)*2*pi/Ndots;

figure();
plot(aa0(k,j)*cos(crcl)-aa0(q,j)*sin(crcl),aa0(k,j)*sin(crcl)+aa0(q,j)*cos(crcl));
hold on;
plot(aa(k,i)*cos(crcl)-aa(q,i)*sin(crcl),aa(k,i)*sin(crcl)+aa(q,i)*cos(crcl),'r');

plot(aa0mf(k,j),aa0mf(q,j),'x');
plot(aamf(k,i),aamf(q,i),'ro');

figure();

plot(aa0mf(k,:),aa0mf(q,:));
hold on;
plot(aamf(k,:),aamf(q,:),'r');


end