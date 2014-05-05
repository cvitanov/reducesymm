function [derival] = zetaderivest(N, np, tp, varargin)

% truncation
derival = 0;
Ap = varargin;
lenAp = length(varargin);

for iter = 1:N
    tmpderival = 0;
    tmpcumsum = cumsum(ones(1,iter))-1;
    idxmax = length(tp)*ones(1,iter) - tmpcumsum(iter:-1:1);  % for each meter index there is a max
    odometer = cumsum(ones(1,iter)); % [1 2 3 4 5 ... iter]
    icount = 0;
    updatedidx = iter;
    odometer_last = odometer;
    updatederiv();
    while 1
        endloop = odometerinc();
        if endloop
            break;
        end
        updatederiv();
    end
    derival = derival + tmpderival * (-1)^iter;
end
derival = -derival;

    function endloop = odometerinc()
        i = iter;
        while i > 0
            if odometer(i) < idxmax(i);
                break;
            end
            i = i-1;
        end
        endloop = 0;
        if i == 0
            endloop = 1;
            return;
        end
        odometer_last = odometer;
        updatedidx = i; %% which index is been updated
        odometer(i:iter) = odometer(i)*ones(1,iter-i+1)+cumsum(ones(1,iter-i+1));
    end
    function updatederiv()
        tplen = sum(np(odometer));
        if  tplen > N
            idxmax(updatedidx) = odometer(updatedidx)-1; %% disregard the combination
            odometer = odometer_last; % roll back odometer
        else
            % odometermat = [odometer;odometermat];
            icount = icount+1;
            % use the term
            tmpval = 1;
            for ai = 1:lenAp
                tmpval = tmpval*sum(Ap{ai}(odometer));
            end
            tmpderival = tmpderival + tmpval*prod(tp(odometer));
        end
    end
end

