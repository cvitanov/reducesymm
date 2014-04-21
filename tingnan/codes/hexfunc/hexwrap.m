function [wx, wy] = hexwrap(x,y,r)

% compute the hex grid index 
[xt, yt] = afftransx(x, y, r);
xhex = floor((floor(xt)+floor(yt)+2)/3);

[xt, yt] = afftransy(x, y, r);
yhex = floor((floor(xt)+floor(yt)+2)/3);

% the basis of hexgrid in cartesian coordinate
v1 = [sqrt(3)*r;0];
v2 = [sqrt(3)*r/2;1.5*r];

arrayv = v1*xhex + v2*yhex;
wx = arrayv(1,:);
wy = arrayv(2,:);
% plot(x, y,'.')
% plot(arrayv(1, :), arrayv(2, :),'o')

function [xt, yt] = afftransx(x, y, r)
AA = [sqrt(3)/3/r,-1/r;
      2*sqrt(3)/3/r,0];
tmp = AA*[x; y];
xt = tmp(1, :);
yt = tmp(2, :);

function [xt, yt] = afftransy(x, y, r)
AA = [ sqrt(3)/3/r, 1/r;
      -sqrt(3)/3/r, 1/r];
tmp = AA*[x; y];
xt = tmp(1, :);
yt = tmp(2, :);
