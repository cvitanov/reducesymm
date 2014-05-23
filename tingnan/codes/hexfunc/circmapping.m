function [ thpp ] = circmapping(r, R1, R2, tph)
%circmapping 
thpp = [];
th = tph(1);
ph = tph(2);

% r = 1;
% w = 0.3;
% 
% R1 = [0;0];
% R2 = [2*r+w;0];
lvec = R2-(R1+r*[cos(th);sin(th)]);
ihat = [cos(th+ph);sin(th+ph)];
proj_lvec = dot(lvec,ihat)*ihat;
npvec = proj_lvec-lvec;
npvec_len = sqrt(dot(npvec,npvec));
if npvec_len > r
    return;
end

nvec = npvec-ihat*sqrt(r*r-dot(npvec,npvec)); 
nhat = nvec/r; % now we get theta
ohat = ihat-2*dot(ihat,nhat)*nhat; % and we get phi

thp = acos(nhat(1));
if nhat(2) < 0
    thp = 2*pi - thp;
end
tmp = cross([nhat;0],[ohat;0]);
php = asin(tmp(3));

thpp = [thp;php];

end

