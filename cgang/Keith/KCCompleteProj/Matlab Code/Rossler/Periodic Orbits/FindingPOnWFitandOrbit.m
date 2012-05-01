% Written by K. Carroll 5/1/2012
load returnmap.txt;
load reformatedxpos.txt
load coefficientsm.txt
N=7;
c = coefficientsm();
sn = returnmap();
xpos = reformatedxpos();
s = size(xpos);
s2 = size(sn);
length = s(1:1,1);
colnum = s(1:1,2);
colnum2 = s2(1:1,2);
snmax = 0;
for i =1:length
    if sn(i:i,1)>snmax
        snmax = sn(i:i,1);
    end
end
snmax;
for j = 1:length
    if sn(j:j,2)==snmax
        scritical = sn(j:j,1);
    end
end
scritical;
numofit = 10000;
sn2 = zeros(numofit,2);
sn4 = zeros(numofit,2);

ds = snmax/(numofit-1);
for j=1:numofit
    sn4(j:j,1) = ds*(j-1);
    sn4(j:j,2) = snp1(ds*(j-1), c);
end

for i=1:numofit
    sn2(i:i,1) = ds*(i-1);
    nextsp = snp1(sn2(i:i,1), c);
    for j=1:N-1
        nextsp = snp1(nextsp, c);
    end
    sn2(i:i,2) = nextsp;
end

sn3 = zeros(length,2);
sn3 = sn2;
sn2(:,2:2) = sn2(:,2)-sn2(:,1);
%find the 0s
rv = 1;
for i=2:numofit
    if sn2(i-1:i-1,2)*sn2(i:i,2)<0
        m = (sn2(i-1:i-1,2)-sn2(i:i,2))/(sn2(i-1:i-1,1)-sn2(i:i,1));
        scross = sn2(i-1:i-1,1)+sn2(i-1:i-1,2)/m;
        sc(rv:rv,1) = scross;
        rv = rv+1;
    end
end
% sc
length2 = size(sc);
PO = zeros(length2, N);
for j=1:length2
    sit = sc(j:j,1);
    for i=1:N
        if sit<scritical
            PO(j:j,i) = 0;
        else
            PO(j:j,i) = 1;
        end
        sit = snp1(sit,c);
    end
end
R = zeros(length2);
for i=1:length2
   v = PO(i,:);
   for j=1:length2
       if (j~=i)&&(R(i:i,1)~=1)
           v2 = PO(j,:);
           a = testifsame(v, v2);
           if a
               R(j:j,1) = 1;
           end
       end
   end 
end
sum = 0;
for i=1:length2
    if R(i:i,1)==1
        sum=sum+1;
    end
end

RedPO = zeros(length2-sum);
rv=1;
for i=1:length2
    if R(i:i,1)~=1
        RedPO(rv:rv,1:N) = PO(i:i,1:N);
        for j=2:length
            if (sn(j-1:j-1,1)<sc(i:i,1))&&(sn(j:j,1)>sc(i:i,1))
                dsn = sn(j:j,1)-sn(j-1:j-1,1);
                t = (sc(i:i,1)-sn(j-1:j-1,1))/dsn;
                vadd = xpos(j-1:j-1,1:colnum)*(1-t)+t*xpos(j:j,1:colnum);
            end
        end
        xn(rv:rv,1:colnum) = vadd;
        rv=rv+1;
    end
end
sn6 = zeros(2*N+1,N);
posit = 1
sn6(1:1,1:2) = [sc(posit:posit,1), sc(posit:1,1)];
sn6(2:2,1:2) = [sc(posit:posit,1), snp1(sc(posit:posit,1),c)];
for i=2:N
    sn6(2*(i-1)+1:2*(i-1)+1,1:2) = [sn6(2*(i-2)+2:2*(i-2)+2,2:2), sn6(2*(i-2)+2:2*(i-2)+2,2:2)];
    sn6(2*(i-1)+2:2*(i-1)+2,1:2) = [sn6(2*(i-1)+1:2*(i-1)+1,1:1), snp1(sn6(2*(i-1)+1:2*(i-1)+1,1:1),c)];
end
sn6(2*N+1:2*N+1,1:2) = sn6(1:1,1:2);
figure(1)
hold all
plot(sn4(:,1), sn4(:,2), 'k', 'linewidth', 2)
% plot(sn4(:,1), sn4(:,1), 'r', 'linewidth', 2)
% plot([scritical, scritical], [snp1(scritical, c), 0], 'k', 'linestyle', '- -', 'linewidth', 2);
plot(sn6(:,1), sn6(:,2),'b', 'linewidth', 2)
xlabel('s_n');
ylabel('s_{n+1}');
dlmwrite('POorder3.txt', PO);
dlmwrite('ReducedPO.txt', RedPO);
dlmwrite('StartingPoints.txt', xn, 'precision', 12);




