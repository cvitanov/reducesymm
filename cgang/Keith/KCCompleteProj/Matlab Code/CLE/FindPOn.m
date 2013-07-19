%Written by K. Carroll 5/1/2012
load reformatedxpos.txt
load returnmap.txt
n = 1;
sn = returnmap();
xpos = reformatedxpos();
s = size(xpos);
length = s(1:1,1);
colnum = s(1:1,2);
sn2 = zeros(length, 2);
numofnonpoints = 0;
for i=1:length
    if sn(i:i,3)==0
        numofnonpoints = numofnonpoints+1;
    end
end
spots = zeros(numofnonpoints,1);
nm = 1;
for i=1:length
    if sn(i:i,3)==0
        spots(nm:nm,1) = i;
        nm = nm+1;
    end
end
length2 = length-numofnonpoints;
sn2 = zeros(length2, 3);
nm = 1;
for i=1:length
%     subv =0;
%     for k=1:numofnonpoints
%         if i>spots(k:k,1)
%             subv = subv+1;
%         end
%     end
    if sn(i:i,3)~=0
        sn2(nm:nm,1:2) = sn(i:i,1:2);
        subv =0;
        for k=1:numofnonpoints
            if sn(i:i,3)>spots(k:k,1)
                subv = subv+1;
            end
        end        
        sn2(nm:nm,3) = sn(i:i,3)-subv;
        for j=1:numofnonpoints
            if sn(i:i,3)==spots(j:j,1)
                sn2(nm:nm,3) = length2+j;
            end
        end
        nm = nm+1;
    end
end
%Now sn2, will have sn, sn+1, then where it goes, if it doesn't have a
%position to go to, it is 0; my goal from here on out is to define where
%these positions should go.  To do this, I will define a matrix which goes
%interpolates where it should be from the given data at each step...In
%other words the column will be the nth iterate, the row will be which
%position it is.

% First i must find a way to find the next sn
places = zeros(numofnonpoints,n-1);
for j=1:length2
    if sn2(j:j,3)>length2
        posp = sn2(j:j,3)-length2;
        nextit = sn2(j:j,2);
        for i=1:n-1
            for k=2:length2
                if (sn2(k-1:k-1,1)<nextit)&&(sn2(k:k,1)>nextit)
                    t = (nextit-sn2(k-1:k-1,1))/(sn2(k:k,1)-sn2(k-1:k-1,1));
                    nextit2 = sn2(k-1:k-1,2)*(1-t)+t*sn2(k:k,2);
                    places(posp:posp, i) = nextit2;
                    nextit = nextit2;
                    break;
                end
            end
        end
    end
end
sn3 = zeros(length2,2);
sn3(:,1:2) = sn2(:,1:2);
for i=1:length2
    sp = sn2(i:i,3);
    for j=1:n-1
        if sp<=length2
            sn3(i:i,2) = sn2(sp:sp,2);
            sp = sn2(sp:sp,3);
        else
            posinpl = sp-length2;
            cnum = n-j;
            sn3(i:i,2) = places(posinpl:posinpl, cnum);
            break;
        end
    end
end
size(places)
sn4 = zeros(length2,2);

sn4(:,1) = sn3(:,1);
sn4(:,2) = sn3(:,2)-sn3(:,1);
rv2 = 0;
for i=2:length2
    if sn4(i-1:i-1,2)*sn4(i:i,2)<=0
        rv2 = rv2+1;
    end
end
xn = zeros(rv2,1:colnum-2);
rv = 1;
for i=2:length2
    if sn4(i-1:i-1,2)*sn4(i:i,2)<=0
        m = (sn4(i-1:i-1,2)-sn4(i:i,2))/(sn4(i-1:i-1,1)-sn4(i:i,1));
        sr = sn4(i-1:i-1,1)-sn4(i-1:i-1,2)/m;
        for j=2:length
            if (sn(j-1:j-1,1)<sr)&&(sn(j:j,1)>sr)
                t = (sr-sn(j-1:j-1,1))/(sn(j:j,1)-sn(j-1:j-1,1));
                posx = j-1;
            end
        end
        xn(rv:rv,1:colnum-2) = xpos(posx:posx,1:colnum-2)*(1-t)+t*xpos(posx+1:posx+1,1:colnum-2);
        rv = rv+1;
    end
end

hold all
plot(sn2(:,1),sn2(:,2));
plot(sn3(:,1), sn3(:,2));
plot(sn4(:,1), sn4(:,2));
dlmwrite('1stiterates.txt', xn,'precision', 12);

