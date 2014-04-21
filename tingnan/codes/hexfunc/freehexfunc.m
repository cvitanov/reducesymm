function [new_tanth, new_tanph, grot] = freehexfunc( th, ph )
% The natural varialbe here are tanth = tan(th) and sin(phi)

%% first region
% (1-sqrt(3)*tan(th))/2 = -(1-sqrt(3)*tan(th))/(1-sqrt(3)*tan(phi))
% tan(phi) = tan(phi-2*pi/3) = (tan(phi)-sqrt(3))/(1+sqrt(3)*tan(phi))
tanphi = tan(ph);
tanth =  tan(th);
if 0.5*sqrt(3)*(tanth + tanphi) >= 1
    
    new_tanth  = (2*(1-sqrt(3)*tanth)./(1-sqrt(3)*tanphi)+1)/sqrt(3);
    new_tanph = (tanphi+sqrt(3))./(1-sqrt(3)*tanphi);
    grot = 2;
    
end

if 0.5*sqrt(3)*(tanth + tanphi) <= -1
    
    new_tanth  =-(2*(1+sqrt(3)*tanth)./(1+sqrt(3)*tanphi)+1)/sqrt(3);
    new_tanph = (tanphi-sqrt(3))./(1+sqrt(3)*tanphi);
    grot = -2;
end


%% second region
% (1-sqrt(3)*tan(th))/2 = (2-sqrt(3)*(tanphi+tanth))/(1+sqrt(3)*tanphi)


if 0.5*sqrt(3)*(tanth + tanphi) < 1 && sqrt(3)*(tanth/2+tanphi) >= 0.5
    
    new_tanth = (1-2*(2-sqrt(3)*(tanphi+tanth))/(1+sqrt(3)*tanphi))/sqrt(3);
    new_tanph = (tanphi-sqrt(3))./(1+sqrt(3)*tanphi);
    grot = 1;
end

if 0.5*sqrt(3)*(tanth + tanphi) > -1 && sqrt(3)*(tanth/2+tanphi) <= -0.5
    
    new_tanth = (2*(2+sqrt(3)*(tanphi+tanth))/(1-sqrt(3)*tanphi)-1)/sqrt(3);
    new_tanph = (tanphi+sqrt(3))./(1-sqrt(3)*tanphi);
    grot = -1;
end

%% 3rd region

if abs(sqrt(3)*(tanth/2+tanphi)) < 0.5
    new_tanth = tanth+2*tanphi;
    new_tanph = tanphi;
    grot = 0;
end


l = (1 - sqrt(3)*new_tanth) / 2;
%%

end

