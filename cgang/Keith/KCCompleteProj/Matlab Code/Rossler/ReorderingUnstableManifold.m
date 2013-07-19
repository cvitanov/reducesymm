%Written by K. Carroll 5/1/2012
load UnstableManifold.txt
xpos = UnstableManifold();

s = size(xpos);
length = s(1:1,1);
colnum = s(1:1,2);
reorder = zeros(length, colnum);

a = 0.2;
b = 0.2;
c = 5.7;
top = (c-sqrt(c^2-4*a*b))/2;
x0n = [top, -(top/a),(top/a)];
trajnum = xpos(1:1, colnum-1);
for i=1:length
    if xpos(i:i,colnum-1)>trajnum
        trajnum = xpos(i:i,colnum-1);
    end
end
initdis = zeros(trajnum, 1);
for i=1:trajnum
    v = xpos(1:1,colnum-2)-x0n;
    dis = dot(v,v);
    initdis(i:i,1) = dis;
end
min = initdis(1:1,1);
spot = 1;
for i=1:trajnum
    min = initdis(1:1,1);
    spot = 1;
    for j=1:trajnum
        if initdis(j:j,1)<min
            spot = j;
            min = initdis(j:j,1);
        end
    end
    reorder(i:i,1:colnum) = xpos(spot:spot,1:colnum);
    xpos(spot:spot,1:colnum) = 10000*ones(1,colnum);
    initdis(spot:spot,1) = 10000;
end
% Now that I have a starting place I can start to find where to put the
% remaining spots
initdis2 = 10000*ones(length, 1);
for i=trajnum+1:length
    if mod(i,100)==0
        i
    end
    for j=1:i-1
        v = xpos(i:i,1:colnum-2)-reorder(j:j,1:colnum-2);
        dis = dot(v,v);
        initdis2(j:j,1) = dis;
    end
    min = initdis2(1:1,1);
    spot = 1;
    for j=1:i-1
        if initdis2(j:j,1)<min
            min = initdis2(j:j,1);
            spot = j;
        end
    end
    reorder = ReplacingVector(reorder, xpos(i:i,1:colnum), spot, i-1);
end

dlmwrite('firstreorder.txt', reorder,'precision', 12);