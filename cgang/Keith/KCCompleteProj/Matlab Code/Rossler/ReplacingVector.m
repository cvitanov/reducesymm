%Written by K. Carroll 5/1/2012
function [xnew] = ReplacingVector(xs, xv, place, endplace)
% reorder = ReplacingVector(reorder, xpos(j:j,1:colnum), spot, i-1);
    sum1 = 0;
    sum2 = 0;
    s = size(xs);
    length=s(1:1,1);
    colnum = s(1:1,2);
    xnew = zeros(length, colnum);
    xrep = zeros(endplace+1, colnum);
    if place~=1&&place~=endplace
        v = xs(place-1:place-1, 1:colnum-2)-xv(1:1,1:colnum-2);
        dis = dot(v,v);
        sum1 = sum1+dis;
        v = xv(1:1,1:colnum-2)-xs(place:place,1:colnum-2);
        dis = dot(v,v);
        sum1 = sum1+dis;
        v = xs(place:place,1:colnum-2)-xs(place+1:place+1,1:colnum-2);
        dis = dot(v,v);
        sum1 = sum1+dis;
        v = xs(place-1:place-1,1:colnum-2)-xs(place:place,1:colnum-2);
        dis = dot(v,v);
        sum2 = sum2+dis;
        v = xs(place:place,1:colnum-2)-xv(1:1,1:colnum-2);
        dis = dot(v,v);
        sum2 = sum2+dis;
        v = xv(1:1,1:colnum-2)-xs(place+1:place+1,1:colnum-2);
        dis = dot(v,v);
        sum2 = sum2+dis;
%         If sum1 is less than sum2, place the new vector between place-1
%         and place; If sum2 is less than sum1, place the new vector between place and place+1 
        if sum1<sum2
            xrep(1:place-1,:) = xs(1:place-1,:);
            xrep(place:place,:) = xv(1:1,:);
            xrep(place+1:endplace+1,:) = xs(place:endplace,:);
        else
            xrep(1:place,:) = xs(1:place,:);
            xrep(place+1:place+1,:) = xv(1:1,:);
            xrep(place+2:endplace+1,:) = xs(place+1:endplace,:);            
        end
    elseif place==1
        xrep(1:1,:) = xs(1:1,:);
        xrep(2:2,:) = xv(1:1,:);
        xrep(3:endplace+1,:) = xs(2:endplace,:);
    elseif place==endplace
        v = xs(endplace:endplace,1:colnum-2)-xs(endplace-1:endplace-1,1:colnum-2);
        d1 = dot(v,v);
        v = xs(endplace-1:endplace-1,1:colnum-2)-xv(1:1,1:colnum-2);
        d2 = dot(v,v);
        if d1<d2
            xrep(1:endplace,:) = xs(1:endplace,:);
            xrep(endplace+1:endplace+1,:) = xv(1:1,:);
        else
            xrep(1:endplace-1,:) = xs(1:endplace-1,:);
            xrep(endplace:endplace,:) = xv(1:1,:);
            xrep(endplace+1:endplace+1,:) = xs(endplace:endplace,:);
        end
    end
    xnew(1:endplace+1, :) = xrep(1:endplace+1,:);
end