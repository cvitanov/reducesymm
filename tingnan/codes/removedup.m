function [ sbmatnew, thmatnew ] = removedup( np, sbmat, thmat )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



%% for period 4 sequences, doing reduction and eliminate all period 1 and period 2 orbits
if np == 4

sbmatnew = [];
thmatnew = [];
for ii = 1:length(sbmat)
    if sbmat(ii, 1) == sbmat(ii, 3) && sbmat(ii, 2) == sbmat(ii, 4)
        continue;
    end
    sbmatnew = [sbmatnew;sbmat(ii,:)];
    thmatnew = [thmatnew;thmat(ii,:)];
end
    return;
end

%% period 6
if np == 6
sbmatnew = [];
thmatnew = [];
for ii = 1:length(sbmat)
    if sbmat(ii, 1) == sbmat(ii, 3) == sbmat(ii, 5) && sbmat(ii, 2) == sbmat(ii, 4) == sbmat(ii, 6)
        continue;
    end
    if sbmat(ii, 1) == sbmat(ii, 4) && sbmat(ii, 2) == sbmat(ii, 5) && sbmat(ii, 3) == sbmat(ii, 6)
        continue;
    end
    sbmatnew = [sbmatnew;sbmat(ii,:)];
    thmatnew = [thmatnew;thmat(ii,:)];
end
return;
end

%% period 8
if np == 8
sbmatnew = [];
thmatnew = [];
for ii = 1:length(sbmat)
    if sbmat(ii, 1) == sbmat(ii, 3) == sbmat(ii, 5) == sbmat(ii, 7) && sbmat(ii, 2) == sbmat(ii, 4) == sbmat(ii, 6) == sbmat(ii, 8)
        continue;
    end
    if sbmat(ii, 1) == sbmat(ii, 5) && sbmat(ii, 2) == sbmat(ii, 6) && sbmat(ii, 3) == sbmat(ii, 7) && sbmat(ii, 4) == sbmat(ii, 8)
        continue;
    end
    sbmatnew = [sbmatnew;sbmat(ii,:)];
    thmatnew = [thmatnew;thmat(ii,:)];
end
return;
end
sbmatnew = sbmat;
thmatnew = thmat;


