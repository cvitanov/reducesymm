
w = 0.3;
r = 1;

Rh = zeros(2,11);

% nearest neighbor
for jj = 0:2:10
    ang = jj*pi/6;
    Rh(:, jj+1) = (2*r+w)*[cos(ang);sin(ang)];
end

% next nearest neight
for jj = 1:2:11
    ang = (jj-1)*pi/6 + pi/6;
    Rh(:, jj+1) = sqrt(3)*(2*r+w)*[cos(ang);sin(ang)];
end


% %% this is the part the generates the permutation of symbols
% allsymbols = 0:11;
% np = 6;
% tmpsblseq = zeros(1, 2);
% NSYMS = 12;
% for ii = 1:NSYMS^np
%     n1 = ii;
%     for jj = 1:np %% determine the digits of the sequence
%         n2 = rem(n1,NSYMS);
%         n1 = fix(n1/NSYMS);
%         tmpsblseq(ii,jj) = n2;
%     end
% end


%%
clear;
np = 6;
tmpsblseq = genNecklaces(np, 12);

% now let us do the symbol reduction for the hard disk 
sblseq = [];
[nx, ~] = size(tmpsblseq);
for ii = 1:nx
    tp = tmpsblseq(ii, :);
    tp = [tp, tp(1)]; % make this periodic to check the last symbol
    test = 1;
    for jj = 2:np+1
        if tp(jj) == tp(jj-1)
            test = 0;
            break;
        end
        numdiff = abs(tp(jj) - tp(jj - 1));
        if numdiff > 6
            numdiff = 12 - numdiff;
        end
        if mod(tp(jj - 1), 2) % odd, change of number should be at least 3
           if numdiff < 3
               test = 0;
               break;
           end
        else
            % even, change of number should be at least 2
           if numdiff < 2
               test = 0;
               break;
           end
        end
    end
    if test == 0
        disp('not a valid sequence for current prunning')
    else
        tp = tp(1:np);
        [~,index] = min(tp(1:np));
        tp = [tp(index:end), tp(1:index-1)];
        sblseq = [sblseq;tp(1:np)];
    end
end

sblseq = unique(sblseq, 'rows');



[sbmat, thmat] = linkminsearch(sblseq);


%%
num = 29
tmpseq = sbmat(num, :);
newth = thmat(num, :);

clear Rv rv
linknum = tmpseq(1)+1;
Rv(:,1) = [0;0];
rv(:,1) = r*[cos(newth(1));sin(newth(1))];


for ii = 2: length(tmpseq)
    Rv(:,ii) = Rh(:, linknum) + Rv(:,ii - 1);
    rv(:,ii) = Rv(:,ii) + r*[cos(newth(ii));sin(newth(ii))];
    linknum = tmpseq(ii)+1;
end
Rv(:,ii+1) = Rh(:, linknum) + Rv(:,ii);
rv(:,ii+1) = Rv(:,ii+1) + r*[cos(newth(1));sin(newth(1))];
linknum = tmpseq(1)+1;
Rv(:,ii+2) = Rh(:, linknum) + Rv(:,ii+1);
rv(:,ii+2) = Rv(:,ii+2) + r*[cos(newth(2));sin(newth(2))];


viscircles(Rv', r*ones(length(Rv),1));
hold on;
plot(Rv(1, :), Rv(2, :), 'o')
plot(rv(1, :), rv(2, :), '-o')
axis image;
% 


%%
% check for intersection to determine the validity of the path
% given a symbol and two theta
jj = 3
symbol = sbmat(num, jj);
R(:,1) = [0;0]; % due to translational symmetry
R(:,2) = Rh(:,symbol+1);

% assume those angles are are provided
th1 = thmat(num, jj);
sndidx = jj+1;
if jj+1 > np
    sndidx = 1;
end
th2 = thmat(num, sndidx);
%
pt1 = R(:,1) + r*[cos(th1);sin(th1)];
pt2 = R(:,2) + r*[cos(th2);sin(th2)];
seg = pt2 - pt1;
seglen = sqrt(dot(seg,seg));
seghat = seg/seglen;
cr(:, 1) = R(:,1) - pt1;
cr(:, 2) = R(:,2) - pt1;

if mod(symbol, 2)
    R(:,3) = Rh(:,mod(symbol-1,12)+1);
    R(:,4) = Rh(:,mod(symbol+1,12)+1);
else
    R(:,3) = Rh(:,mod(symbol-2,12)+1);
    R(:,4) = Rh(:,mod(symbol+2,12)+1);
end

cr(:, 3) = R(:,3) - pt1;
cr(:, 4) = R(:,4) - pt1;
test = 0;
projm = zeros(2,4);
for ii = 1:4
    proj = dot(cr(:, ii),seghat);
    if proj < 0 || proj > seglen
        continue;
    else
       projv = proj*seghat;
       distv = pt1 + projv - R(:,ii);
       projm(:,ii) = projv+pt1;
       distl = sqrt(dot(distv, distv));
       if distl < r
           test  = 1;
           break;
       end
    end
end

% close gcf
viscircles(R', r*ones(length(R),1));
line([pt1(1),pt2(1)], [pt1(2),pt2(2)]); hold on;
plot(projm(1,:),projm(2,:),'o')
axis image

%% for period 4 sequences, doing reduction and eliminate all period 1 and period 2 orbits

sbmatnew = [];
thmatnew = [];
for ii = 1:length(sbmat)
    if sbmat(ii, 1) == sbmat(ii, 3) && sbmat(ii, 2) == sbmat(ii, 4)
        continue;
    end
    sbmatnew = [sbmatnew;sbmat(ii,:)];
    thmatnew = [thmatnew;thmat(ii,:)];
end
sbmat = sbmatnew;
thmat = thmatnew;

%%

sbmatnew = [];
thmatnew = [];
for ii = 1:length(sbmat)
    if sbmat(ii, 1) == sbmat(ii, 3) && sbmat(ii, 2) == sbmat(ii, 4)
        continue;
    end
    sbmatnew = [sbmatnew;sbmat(ii,:)];
    thmatnew = [thmatnew;thmat(ii,:)];
end
sbmat = sbmatnew;
thmat = thmatnew;

