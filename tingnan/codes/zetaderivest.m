function [derival] = zetaderivest(N, np, tp, varargin)

% truncation
derival = 0;
Ap = varargin;
lenAp = length(varargin);
[C, ia, ~] = unique(np);
npmax = [ia(2:end)-1; length(np)]';
itermax = (N-mod (N, min(C)))/min(C);
for iter = 1:itermax
    tmpderival = 0;
    tmpcumsum = cumsum(ones(1,iter))-1;
    odometer = cumsum(ones(1,iter)); % [1 2 3 4 5 ... iter]
    idxmax = length(tp)*ones(1,iter) - tmpcumsum(iter:-1:1);  % for each meter index there is a max
    icount = 0;
    setodomax();
    updatederiv();
    while 1
        endloop = odometerinc();
        if endloop
            break;
        end
        updatederiv();
    end
    derival = derival + tmpderival * (-1)^iter;
    iter
    tmpderival
    icount
end
derival = -derival;

    function endloop = odometerinc()
        i = iter;
        while i > 0
            remlength = N - sum(np(odometer(1:i-1)));
            if odometer(i) < min(idxmax(i), sum(npmax(C == remlength)));
                break;
            end
            i = i-1;
        end
        endloop = 0;
        if i == 0
            endloop = 1;
            return;
        end
        odometer(i:iter) = odometer(i)*ones(1,iter-i+1)+cumsum(ones(1,iter-i+1));
    end
    function updatederiv()
        tplen = sum(np(odometer));
        if  tplen > N
        else
            % odometermat = [np(odometer);odometermat];
            icount = icount+1;
            % use the term
            tmpval = 1;
            for ai = 1:lenAp
                tmpval = tmpval*sum(Ap{ai}(odometer));
            end
            tmpderival = tmpderival + tmpval*prod(tp(odometer));
        end
    end
    function setodomax()
        divnum = N - iter * min(C);
        longestcycle = divnum + min(C);
        idxmax = npmax(C==longestcycle)*ones(1,iter) - tmpcumsum(iter:-1:1);
    end
end

