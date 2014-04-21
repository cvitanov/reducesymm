function [ xyp ] = periodicorbit(func, bdfunc, xy0, d, n)
%given an initial guess of x0, find the exact position of a periodic orbit
%of length n

toly = 1e-8;
tolx = 1e-8;
tolm = 1e-8;
stepmax = 5e-2;
niter = 50;
dim = d;
%
% construct the function
% xy(1), xy(2), xy(3), xy(4), ... xy(n) = xy(1)
xy = zeros(2, n);
xy(:, 1) = xy0;
if n >=2
    for ii = 2:n
        xy(:, ii) = func(xy(:, ii-1));
    end
end
xvec = reshape(xy, [dim*n 1]);
yvec = ndshootingfunction(xvec);
jj = 1;
while jj <= niter
    [jac, ~] = jacobianest(@ndshootingfunction, xvec);
    [L, U, ~] = lu(jac);
    if abs(prod(diag(U))) < 1e-8
        warning('jacobian is close to singular! stop iteration');
        xyp = [];
        break;
    end
    delta = -jac\yvec;
    gradt = yvec'*jac;
    steplen = sqrt(dot(delta, delta));
    % this is to avoid very long shoot that unfortunately get out of bound immediately
    if steplen > stepmax
        delta = delta * stepmax/steplen;
    end
    [lambda, ynew, fnew, check] = lnsrch(xvec, yvec, delta);
    xvec = xvec+lambda*delta;
    yvec = ynew;
    
    if max(abs(yvec)) < toly
        disp('find a solution');
        xyp = reshape(xvec, [dim, n]);
        break;
    end
    
    if check == 1
        den = max(fnew, 0.5*n);
        test = max(abs(gradt)/den);
        if test < tolm
         disp('check')
            xyp = [];
            break;
        end       
    end
    
    % when no solution is found, check for boundary;
    if bdfunc(xvec)
        xyp = [];
        disp('out of boundary! stop iterations');
        break;
    end
    jj = jj + 1;
end

if jj >= niter
    disp('not converging');
%     figure; hold on;
%     seq
%     plot(seq(1,1), seq(2,1), 'o');
%     plot(seq(1,:), seq(2,:), '.-');
%     plot(seq(1,end), seq(2,end), '+');
%     val = diff(seq, 1, 2);
%     figure;
%     plot(sqrt(val(1, :).^2 + val(2, :).^2));
    xyp = [];
%     error('abort');
end

    function [ out ] = ndshootingfunction(xx)
        y = reshape(xx, [dim n]);
        if n>=2
            for ii = 2:n
                tmp(:, ii) = y(:, ii) - func(y(:, ii-1));
            end
        end
        tmp(:, 1) = y(:, 1) - func(y(:, n));

        out = reshape(tmp, [dim*n 1]);
    end
    
    function [lam, ynew, fnew, check] = lnsrch(xold, yold, detx)
        fold = qualfunc(yold);
        localslope = -2*fold;
        alpha = 1e-4;
        lamin = tolx/steplen; %%
        lam = 1;
        check = 0;
        while 1
            xnew = xold+lam*detx;
            ynew = ndshootingfunction(xnew);
            fnew = qualfunc(ynew);
            if lam < lamin
                check = 1;
                return;
            end
            if fnew < fold+alpha*lam*localslope
                return;
            end
            if lam == 1
                tmplam = -0.5*localslope/(fnew-fold-localslope); % fold/(fold+fnew);
            else
                % fnew has been updated.
                lam1  = lam;
                fnew1 = fnew;

                rhs =  [fnew1-fold-lam1*localslope
                        fnew2-fold-lam2*localslope];
                lhs = [    1/lam1^2,   -1/lam2^2;
                       -lam2/lam1^2, lam1/lam2^2]*rhs/(lam1-lam2);
                a = lhs(1);
                b = lhs(2);
                if abs(a) < 1e-8;
                    tmplam = -localslope/(2*b);
                else
                    disc = b^2-3*a*localslope;
                    if disc < 0
                        tmplam = 0.5*lam1;
                    elseif b < 0
                        tmplam = (-b+sqrt(disc))/(3*a); 
                    else
                        tmplam = -localslope/(b+sqrt(disc));
                    end
                end
                tmplam = min(0.5*lam1, tmplam);
            end
            lam2 = lam;
            fnew2 = fnew;
            lam = max(tmplam, 0.1*lam);
        end
    end

    function [out] = qualfunc(fxy)
        out = 0.5*dot(fxy,fxy);
    end

end

