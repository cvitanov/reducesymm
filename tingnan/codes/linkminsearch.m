function [ sbmat, thmat ] = linkminsearch( sblseq, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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

[nx, ~] = size(sblseq);
sbmat = [];
thmat = [];
for kk = 1:nx
    tmpseq = sblseq(kk, :);    
    th = rand(1, length(tmpseq))*2*pi;
    if(~isempty(varargin))
        th = varargin{1}(kk, :);
    end
    exitflag = 0;
    [newth, ~, exitflag] = fminsearch(@linkfunc, th,  optimset('TolX',1e-14));
    % while(~exitflag)
    %     [newth, ~, exitflag] = fminsearch(@linkfunc, th,  optimset('TolX',1e-10));
    % end
    newth = mod(newth,2*pi);
    thseq = [newth newth(1)];
    for jj = 1:length(newth)
        iftest = testlink(sblseq(kk, jj), thseq(jj), thseq(jj+1));
        if iftest == 1
            break;
        end
    end
    if iftest == 0
        % accept this as a solution
        sbmat = [sbmat;sblseq(kk,:)];
        thmat = [thmat;newth];
    end
end


    function [test] = testlink(symbol, th1, th2)
        R(:,1) = [0;0]; % due to translational symmetry
        R(:,2) = Rh(:,symbol+1);
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
%         projm = zeros(2,4);
        for ii = 1:4
            proj = dot(cr(:, ii),seghat);
            if proj < 0 || proj > seglen
                continue;
            else
                projv = proj*seghat;
                distv = pt1 + projv - R(:,ii);
%                 projm(:,ii) = projv+pt1;
                distl = sqrt(dot(distv, distv));
                if distl < r
                    test  = 1;
                end
            end
        end
%         if test == 0
%             close gcf
%             viscircles(R', r*ones(length(R),1));
%             line([pt1(1),pt2(1)], [pt1(2),pt2(2)]); hold on;
%             plot(projm(1,:),projm(2,:),'o')
%             axis image
%             waitforbuttonpress;
%         end
    end

    function [linklen] = linkfunc(thvec)
        linklen = 0;
        thvectmp = [thvec,thvec(1)];
        for ii = 1:length(thvec)
            sbidx = tmpseq(ii)+1;
            dispv = [cos(thvectmp(ii+1));sin(thvectmp(ii+1))] - [cos(thvectmp(ii));sin(thvectmp(ii))] + Rh(:, sbidx);
            linklen = linklen + sqrt(dot(dispv, dispv));
        end
    end
       
end

