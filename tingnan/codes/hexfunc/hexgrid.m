function hexgrid(r, rr)

nx = 3;
ny = 3;

v=30:60:390;
cv=r*cosd(v);
sv=r*sind(v);

for m=-ny:ny
    for k=-2*nx:2*nx
        xseries = k*sqrt(3)*r+cv;
        yseries = m*3*r+sv;
        viscircles([k*sqrt(3)*r,m*3*r],rr, 'LineWidth', .5, 'LineStyle', '--', 'EdgeColor', 'b');
%         [xt, yt] = afftransy(xseries, yseries, r);
        line(xseries,yseries,'tag','h');
    end
end

for m=-ny:ny
    for k=-2*nx:2*nx
        xseries = sqrt(3)*r/2 + k*sqrt(3)*r+cv;
        yseries = 1.5*r+m*3*r+sv;
        viscircles([sqrt(3)*r/2+k*sqrt(3)*r,1.5*r+m*3*r],rr, 'LineWidth', .5, 'LineStyle', '--', 'EdgeColor', 'b');
%         [xt, yt] = afftransy(xseries, yseries, r);
        line(xseries,yseries,'tag','h');
    end
end
hold on;
axis equal

% x = 10*r*rand(1,10) - 5*r;
% y = 10*r*rand(1,10) - 5*r;
% 
% [wx wy] = hexwrap(x,y,r);
% % [xt, yt] = afftransx(x, y, r);
% % xhex = floor((floor(xt)+floor(yt)+2)/3);
% % 
% % [xt, yt] = afftransy(x, y, r);
% % yhex = floor((floor(xt)+floor(yt)+2)/3);
% % 
% % xhex
% % yhex
% % v1 = [sqrt(3)*r;0];
% % v2 = [sqrt(3)*r/2;1.5*r];
% % 
% % arrayv = v1*xhex + v2*yhex;
% % arrayv
% plot(x, y,'.')
% plot(wx, wy,'o')

% function [xt, yt] = afftransx(x, y, r)
% AA = [sqrt(3)/3/r,-1/r;
%       2*sqrt(3)/3/r,0];
% % AA = [1,0;
% %       0,1];
% tmp = AA*[x; y];
% xt = tmp(1, :);
% yt = tmp(2, :);
% 
% function [xt, yt] = afftransy(x, y, r)
% AA = [ sqrt(3)/3/r, 1/r;
%       -sqrt(3)/3/r, 1/r];
% % AA = [1,0;
% %       0,1];
% tmp = AA*[x; y];
% xt = tmp(1, :);
% yt = tmp(2, :);
