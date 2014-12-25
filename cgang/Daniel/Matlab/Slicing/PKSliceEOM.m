
% Equations of motion for the Porter Knobloch system on a slice
% Takes a system state y = {x1,x2,y1,y2,z} and returns the velocity
% vector dy = [vhat phidot] as per eqs. 10.49 and 10.50 in Chaosbook
function dy = PKSliceEOM(t,y,params)

mu1 = params{1};
mu2 = params{2};
a1 = params{3};
a2 = params{4};
b1 = params{5};
b2 = params{6};
c1 = params{7};
c2 = params{8};
e2 = params{9};

if params{12} == 1
    xprime = params{10}; % Set template 1 as slicing template
elseif params{12} == 2
    xprime = params{11}; % Set template 2 as slicing template
end

T = TPK;

tgprime = T*xprime(:); % Calculate group tangent for template as a column vector

x1 = y(1);
x2 = y(2);
y1 = y(3);
y2 = y(4);

x = [x1 x2 y1 y2]';

tg = T*x; % Calculate group tangent for current state of PK system

v = zeros(4,1); 

% Porter Knobloch equations in the full Cartesian space
v(1) = a1*x1^3 + b1*x1*x2^2 + c1*x1*x2 + a1*x1*y1^2 + b1*x1*y2^2 + mu1*x1 + c1*y1*y2;
v(2) = a1*x1^2*y1 + c1*x1*y2 + b1*x2^2*y1 - c1*x2*y1 + a1*y1^3 + b1*y1*y2^2 + mu1*y1;
v(3) = c2*(x1^2 - y1^2) + e2*y2 + mu2*x2 + x2*(a2*x1^2 + a2*y1^2) + 2*b2*x2*y2^2 + b2*x2*(x2^2 - y2^2);
v(4) = a2*x1^2*y2 + 2*c2*x1*y1 + b2*x2^2*y2 - e2*x2 + a2*y1^2*y2 + b2*y2^3 + mu2*y2;
 
phidot =  ((v.')*tgprime)/(tg.'*tgprime); % Calculate phase velocity by eq. 10.49 in Chaosbook

vhat = v - phidot*tg; % Calculate velocity perpendicular to the slice by eq. 10.50 in Chaosbook

dy = [vhat(1), vhat(2), vhat(3), vhat(4), phidot].';