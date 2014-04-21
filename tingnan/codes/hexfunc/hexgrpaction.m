function [th_p, ph_p, grot] = hexgrpaction(th, ph)
% group action
% grot, a number that referes to how many clockwise pi/3 rotation is
% performed. 

valrem = mod(th, pi/3);
valdiv = floor(th/pi*3);

th_p = valrem;
ph_p = mod(ph - valdiv * pi/3, 2*pi); % rotate first

if ph_p > pi
    ph_p = ph_p - 2*pi;
end


grot = valdiv;
% th
% grot

if th_p > pi/6
    % back rotation
    th_p = th_p - pi/3;
    ph_p = ph_p - pi/3;
    grot = grot + 1;
end

% (th - th_p)*180/pi
% (ph - ph_p)*180/pi
end

