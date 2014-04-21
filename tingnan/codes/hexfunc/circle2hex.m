function [th_p, ph_p, grot] = circle2hex(th, ph, r)
% when leaving the scattering center, we now need to compute the point on
% the Hex boundary


%% bugs find here. ll
if abs(ph) <= pi/2
ll = r*sin(th) + (sqrt(3)/2 - r*cos(th))*tan(ph); 
    if abs(ll) <= 0.5
        tanthp = (2*ll)/sqrt(3);
        th_p = atan(tanthp);
        ph_p = ph;
        grot = 0;
        return;
    end
end


if ph >= 0    
    ph_0 = atan((0.5-r*sin(th)) / (sqrt(3)/2 - r*cos(th)));
    l1 = (sqrt(3)/2-r*cos(th)) / cos(ph_0);
    l2 = l1*sin(ph-ph_0) / cos(ph - pi/3);
    tanthp = (2*l2-1)/sqrt(3);
    th_p = atan(tanthp);
    ph_p = ph - pi/3;
    grot = 1;
end

if ph < 0
    ph_0 = atan((0.5+r*sin(th)) / (sqrt(3)/2 - r*cos(th)));
    l1 = (sqrt(3)/2-r*cos(th)) / cos(ph_0);
    l2 = -l1*sin(ph+ph_0) / cos(ph + pi/3);
    tanthp = (1-2*l2)/sqrt(3);
    th_p = atan(tanthp);
    ph_p = ph + pi/3;    
    grot = -1;
end
% 
% if abs(th_p) > pi/6
% 
%     th * 180/pi
%     ph * 180/pi
%     ph_0 * 180/pi
%     l2
%     th_p * 180/pi
%     ph_p * 180/pi
% end

end


