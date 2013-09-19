function [value,isterminal,direction] = PKevents(t,x,params)

T = TPK;
xp1 = params{10}; % Set template 1
xp2 = params{11}; % Set template 2

if params{12} == 1
    xprime = xp1; % Set template 1 as slicing template
elseif params{12} == 2
    xprime = xp2; % Set template 2 as slicing template  
end

t1 = T*xp1(:); % Calculate group tangent for template 1 as a column vector
t2 = T*xp2(:); % Calculate group tangent for template 1 as a column vector
tp = T*xprime(:); % Calculate group tangent for template currently in use
X = x(1:4);
tg = T*X(:); % Calculate group tangent for current point

tg = tg/norm(tg);
tp = tp/norm(tp);

% Define stop condition when <(x-template)|t> switches signs
value = [dot(x(1:4),t1(:));...
         dot(x(1:4),t2(:));...
         dot(tg(:),tp(:))]; % value(3) checks that group tangents for point and current template are not orthogonal

if params{12} == 1
    value(1) = 0;
end
if params{12} == 2
    value(2) = 0;
end


isterminal = [1; 1; 1];  % Force integrator to stop

direction = [1; 1; 0]; % Stop both when <(x-template)|t> goes from - to + and from + to -

end % end events 