function [xyp, xyt] = relaxation2d( func, bdfunc, xy0, n )
% The relaxation method for finding periodic orbits

C{1} = [ 1  0; 
       0  1];
C{2} = [-1  0; 
       0  1];
C{3} = -C{2};

% C{4} = [ 0 -1; 
%       -1  0];
% C{5} = -C{4};

tau = 1e-3;
tol = 1e-5;
bnd = 2;
S = cell(3, 1);
cmaps = hsv(3);
nn = 1e4;
xyp = [];
xyt = [];
for jj = 1:3
    S{jj} = @(xypt) xypt + tau*C{jj}*(ifunc(xypt) - xypt);
    xy = zeros(2, 1);
    xy(:, 1) = xy0;
    ii = 2;
    converge = 0;
    while 1
        for kk = ii:ii+10
            xy(:, kk) = S{jj}(xy(:, kk-1));
        end
        
        if bdfunc(xy(:, end)) || ii > nn
            converge = 0;
            break;
        end

        if pdist(xy(:, ii-1:ii)') < tol 
            converge = 1;
            break;
        end
        
        ii = kk + 1;
    end
    if converge == 1
%         plot(xy(1,1), xy(2,1), 'o', 'color', cmaps(jj, :));
%         plot(xy(1,:), xy(2,:), '.-', 'markersize', 0.5, 'linewidth', 0.5);
%         plot(xy(1,end), xy(2,end), '+', 'color', cmaps(jj, :));
        newbdfunc = @(xypt) bdfunc(xypt) || pdist([xypt(1),xypt(2);xy(1,end),xy(2,end)]) > 1e-1;
        tmpxyp = periodicorbit(func, newbdfunc, [xy(1,end); xy(2,end)], 2, n);
        if ~isempty(tmpxyp)
            figure(1);
            hold on;
%             plot(xy(1,:),  xy(2,:), '.', 'markersize', 1, 'color', cmaps(jj, :)); 
%             plot(xy(1,end), xy(2,end), '+', 'color', cmaps(jj, :));
            plot(tmpxyp(1,1), tmpxyp(2,1), 'square', 'color', cmaps(jj, :));
            xyp = [xyp, tmpxyp];
        else
            xyt = [xyt, [xy(1,end); xy(2,end)]];
        end
    end
end

    function xy1 = ifunc(xy0)
        xy1 = func(xy0);
        if n > 1
            for it = 2:n
                xy1 = func(xy1);
            end
        end
    end

end