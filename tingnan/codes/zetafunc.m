function [ zetacoefs ] = zetafunc(np, tp)
%
Np = max(np);
zetacoefs = [1,zeros(1, Np)];
for ii = 1:length(tp)
    maxidx = Np - np(ii) + 1;
    zetacoefs = zetacoefs - tp(ii)*[zeros(1, Np + 1 - maxidx), zetacoefs(1:maxidx)];
end
end

