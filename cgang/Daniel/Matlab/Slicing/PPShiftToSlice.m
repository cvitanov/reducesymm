
function [THETA,fval,xhat] = PPShiftToSlice(theta0,x,template) 
options = optimset('Display','off','TolFun',3e-14);
[THETA,fval] = fsolve(@nestedfun,theta0,options);
% Nested function that computes the objective function     
    function y = nestedfun(theta)
        T = TPK; % Calculate generator of infinitsemal rotations 
        tgprime = T*template(:); % Calculate group tangent for template
        x = x(:); %Make sure x is a column vector
        G = gPK(theta);
        y = (G\x).'*tgprime;     
    end
G = gPK(THETA);
xhat = (G\x);
end