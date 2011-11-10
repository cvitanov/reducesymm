i=59; j=34; % choose snapshots on 2 periodic orbits
k=5;q=6; % projection

Ndots=100;
crcl=(0:Ndots-1)*2*pi/Ndots;

% figure();
plot(aa0(k,j)*cos(crcl)-aa0(q,j)*sin(crcl),aa0(k,j)*sin(crcl)+aa0(q,j)*cos(crcl));
hold on;
plot(aa(k,i)*cos(crcl)-aa(q,i)*sin(crcl),aa(k,i)*sin(crcl)+aa(q,i)*cos(crcl),'r');

plot(aa0mf(k,j),aa0mf(q,j),'x');
plot(aamf(k,i),aamf(q,i),'ro');