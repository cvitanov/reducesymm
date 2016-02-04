function [xy, uv] = LorentzGasChaoticBouncing(x0, v0, r, w)
% The function generates a trajectory in the Lorentz gas system with the
% given initial condition [x0, v0] and the parameter [r w]
restitution = 0.0;
a_mat = [2*r+w,r+w/2;0,(r+w/2)*sqrt(3)];

t_span = 0:1e-3:10;
ode_opts = odeset('Events', @events, 'MaxStep', 1e-2);
xy = x0';
ve = v0;
uv = v0';
for it = 1:50
    [T,Y,TE,YE] = ode45(@odedy,t_span,[xy(end,:)';ve],ode_opts);
    [distances,centers] = ComputeDist(YE);
    [~,index] = min(distances);
    disp = (YE(1:2)'-centers(:,index))/r;
    ve = YE(3:4)';
    vn = disp*dot(disp, ve);
    ve = (ve-(1 + restitution)*vn);ve-2*vn;
    ve = ve / sqrt(dot(ve, ve));
    xy = [xy;YE(1:2)];
    uv = [uv;ve'];
end

line_dist = sqrt(3)*(r+w/2)*2;
for ii = -5:5
    CreateCircles(0,ii*line_dist);
    CreateCircles(r+w/2,(ii+0.5)*line_dist);
end
axis image;
plot(xy(:,1),xy(:,2)); hold on;
quiver(xy(:,1),xy(:,2),uv(:,1),uv(:,2), 0.2);

    function CreateCircles(x0, y0)
        xpts = x0+(2*r+w)*(-10:10);
        ypts = y0+zeros(1,length(xpts));
        for jj = 1:length(xpts)
            rectangle('Position',[xpts(jj)-r,ypts(jj)-r,2*r,2*r], 'Curvature', [1 1], 'FaceColor',[.5,.5,.5], 'EdgeColor', [0,0,0]); hold on;
        end
    end


    function dydt = odedy(t,y)
        dydt = [y(3);y(4);0;0];
    end
    function [value,isterminal,direction] = events(t,y)
        dists = ComputeDist(y);
        value = min(dists)-r;
        isterminal = 1;
        direction = -1;
    end
    function [dists,centers] = ComputeDist(y)
        a1_a2 = floor(linsolve(a_mat, [y(1);y(2)]));
        centers = a_mat*(repmat(a1_a2,[1,4])+[0 0 1 1; 0 1 0 1]);
        disps = repmat([y(1);y(2)],[1,4])-centers;
        dists = sqrt(dot(disps,disps,1));
    end
end

