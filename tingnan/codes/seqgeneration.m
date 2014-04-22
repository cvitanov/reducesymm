function [ symbolseqs ] = seqgeneration(symbolist, np)
symbolist = 0:11;
nsyms = length(symbolist);
tmpsblseq = zeros(1, nsyms);
hasharray = [];
for ii = 1:nsyms^np
    n1 = ii;
    for jj = 1:np %% determine the digits of the sequence
        n2 = rem(n1,nsyms);
        n1 = fix(n1/nsyms);
        tmpsblseq(ii,jj) = symbolist(n2+1);
    end
    % now we created a symbol array
end

 
%% now let us do the symbol reduction:
% sblseq = [];
% for ii = 1:nsyms^np
%     tp = tmpsblseq(ii, :);
%     tp = [tp, tp(1)]; % make this periodic to check the last symbol
%     test = 1;
%     for jj = 2:np+1
%         if tp(jj) == tp(jj-1)
%             test = 0;
%             break;
%         end
%         numdiff = abs(tp(jj) - tp(jj - 1));
%         if numdiff > 6
%             numdiff = 12 - numdiff;
%         end
%         if mod(tp(jj - 1), 2) % odd, change of number should be at least 3
%            if numdiff < 3
%                test = 0;
%                break;
%            end
%         else
%             % even, change of number should be at least 2
%            if numdiff < 2
%                test = 0;
%                break;
%            end
%         end
%     end
%     if test == 0
%         disp('not a valid sequence for current prunning')
%     else
%         tp = tp(1:np);
%         [~,index] = min(tp(1:np));
%         tp = [tp(index:end), tp(1:index-1)];
%         sblseq = [sblseq;tp(1:np)];
%     end
% end

symbolseqs = unique(sblseq, 'rows');

    function [] = canonicalform(sblrow)
        % generate all 12 permutations
        % actually we do not need 12
        minelem = min(sblrow);
        minidxs = (sblrow == minelem);
        
    end

end

