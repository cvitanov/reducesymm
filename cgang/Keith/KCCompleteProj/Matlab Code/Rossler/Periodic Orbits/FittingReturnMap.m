% Written by K. Carroll 5/1/2012
load returnmap.txt;

sn = returnmap();

s = size(sn);
length = s(1:1,1);
colnum = s(1:1,2);
numoftraj = sn(1:1,colnum-1);
numofnonpoints = 0;
for i=1:length
    if sn(i:i,colnum-2)==0
        numofnonpoints = numofnonpoints+1;
    end
    if sn(i:i,colnum-1)>numoftraj
        numoftraj = sn(i:i,colnum-1);
    end
end
length2 = length - numofnonpoints;
sn2 = zeros(length2, 2)
rv = 1;
for i=1:length
    if sn(i:i,colnum-2)~=0
        sn2(rv:rv,1:2) = sn(i:i,1:2);
        rv = rv+1;
    end
end
power = 8;

S = zeros(2*power-1,1);

P = zeros(power,1);
P2 = zeros(power,1);

coeff = zeros(power,1);

comatrix = zeros(power,power);
solution = zeros(power,power);
% length = 3;
for i=1:length2
    for j = 1:(2*power-1)
        S(j:j,1) = S(j:j,1) + sn2(i:i,1)^(j-1);
    end
    for k = 1:power
        P(k:k,1) = P(k:k,1) + sn2(i:i,2)*sn2(i:i,1)^(k-1);
    end
end

for i=1:power
    for j=1:power
        solution(i:i,j) = S(i+j-1:i+j-1,1);
        solution2(i:i,j) = S(i+j-1:i+j-1,1);
    end
end

comatrix = solution;

bot = det(solution);
bot2 =det(solution2);
for i=1:power
    comatrix = solution;
%     comatrix2 = solution2;
    comatrix(:,i) = P;
%     comatrix2(:,i) = P2;

    top = det(comatrix);
%     top2 = det(comatrix2);
    coeff(i:i,1) = top/bot;
%     coeff(i:i,2) = top2/bot2;
end
newdata = zeros(length2,1);
newdata2 = zeros(length2,1);

for i=1:length2
    spot(i:i,1) = i;
    for j=1:power
        newdata(i:i,1) = sn2(i:i,1);
        newdata2(i:i,1) = newdata2(i:i,1) + coeff(j:j,1)*(sn2(i:i,1))^(j-1);
    end
end

hold all;
plot(sn2(:,1),sn2(:,2),'k', 'linestyle', 'none', 'Marker', 'o','MarkerFaceColor','k', 'MarkerSize', 5);
% plot(r,z, 'b', 'linewidth', 4);
plot(newdata, newdata2, 'g','linewidth',3);

dlmwrite('coefficientsm.txt', coeff);