%Written by K. Carroll 5/1/2012
load firstreorder.txt
xpos = firstreorder();

a = 0.2;
b = 0.2;
c = 5.7;
top = (c-sqrt(c^2-4*a*b))/2;
x0n = [top, -(top/a),(top/a)];

s = size(xpos);
length = s(1:1,1);
colnum = s(1:1,2);

numoftraj = xpos(1:1,colnum-1);
for i=1:length
    if xpos(i:i,colnum-1)>numoftraj
        numoftraj = xpos(i:i,colnum-1);
    end
end
crossovers = length/numoftraj;
crsmat = zeros(numoftraj, crossovers);
for i=1:numoftraj
    placing = 1;
    for j=1:length
        if xpos(j:j,colnum-1)==i
            crsmat(i:i,placing:placing) = xpos(j:j,colnum);
            placing = placing+1;
        end
    end
end
turningpoint = zeros(numoftraj,1);
for i=1:numoftraj
    spotprev = crsmat(i:i,1:1);
    for j=1:crossovers
        for k=1:crossovers
            if crsmat(i:i,k:k) == j
                spot = k;
            end
        end
        if spot<spotprev
            break;
        else
            spotprev = spot;
        end
    end
    turningpoint(i:i,1) = crsmat(i:i,spot:spot);
end
% Now turningpoint will tell me when to stop counting each trajectory
subtractvalue = 0;
for i=1:numoftraj
    subtractvalue = subtractvalue+(crossovers-turningpoint(i:i,1));
end
length2 = length-subtractvalue;
xpos2 = zeros(length2, colnum);
% sn2 = zeros(length2,4);
sn1 = zeros(length2,5);
sn3 = zeros(length2,2);
posit = 1;
for i=1:length
    traj = xpos(i:i,colnum-1);
    if xpos(i:i,colnum)<=turningpoint(traj:traj,1)
        xpos2(posit:posit,1:colnum) = xpos(i:i,1:colnum);
        posit = posit+1;
    end
end
% Now I must find sn1 and then sn+1
for i=1:numoftraj
    for j=1:turningpoint(i:i,1)-1
        for k=1:length2
            if (xpos2(k:k,colnum-1)==i)&&(xpos2(k:k,colnum:colnum)==j)
                place = k;
            end
        end
        for k=1:length2           
            if (xpos2(k:k,colnum-1)==i)&&(xpos2(k:k,colnum:colnum)==j+1)
                sn1(place:place,3)= k;
                sn1(place:place,4) = i;
                sn1(place:place,5) = 1;
            end
        end
    end
    i
end
v = xpos2(1:1,1:colnum-2)-x0n(1:1,:);
dis = dot(v,v);
sum = 0;
sum = sum+sqrt(dis);
sn1(1:1,1) = sum;
for i=2:length2
    v = xpos2(i:i,1:colnum-2)-xpos2(i-1:i-1,1:colnum-2);
    dis = dot(v,v);
    sn1(i:i,1) = sn1(i-1:i-1,1)+sqrt(dis);
end
for i=1:length2
    spot = sn1(i:i,3);
    sn3(i:i,1:2) = [sn1(i:i,1), sn1(i:i,1)];
    if spot~=0
        sn1(i:i,2) = sn1(spot:spot,1);
    end
end

hold all
plot(sn1(:,1), sn1(:,2),  'r', 'linestyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 3);
plot(sn3(:,1), sn3(:,2), 'b');

dlmwrite('reformatedxpos.txt', xpos2, 'precision', 12)
dlmwrite('returnmap.txt', sn1, 'precision', 12)