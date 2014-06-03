function [ newtable ] = symbolConvert( sbtable )
% convert disk labels to fundamental domain symbols
% from Shreiber
% Symbols	Last long	Last short	Next same	Next opposite
% a	1	2	x	
% b	3	4	x	
% c	5	6	x	
% d	5	4		x
% e	3	2		x
% f	1	-		x
% A	2	1	x	
% B	4	3	x	
% C	6	5	x	
% D	4	5		x
% E	2	3		x
% F	-	1		x


fdtable = [
1   2   1   0	
3	4	1   0	
5	6	1	0
5	4	0	1
3	2	0	1
1	0	0	1
2	1	1	0
4	3	1	0
6	5	1	0
4	5	0	1
2	3	0	1
0	1	0	1];
strlist = 'abcdefABCDEF';
NSYMS = 12;

[nx, ny] = size(sbtable);
newtable = '';
for ii = 1:nx
    currline = [sbtable(ii,:), sbtable(ii,1)];
    currstr = [];
    curridx = [];
    aofchange = zeros(1,ny+1);
    dofchange = zeros(1,ny+2);
    for jj = 2:ny+1 % all the relative direction change
        % if the change is greater than 6, mark as 1, less than 6 mark as 0
        lastchange = mod(currline(jj)-currline(jj-1),12);
        aofchange(jj) = lastchange;
        dofchange(jj) = 0;
        if lastchange > 6
            dofchange(jj) = 1;
            aofchange(jj) = 12 - lastchange;
        end
    end
    dofchange(ny+2) = dofchange(2);
    for jj = 2:ny+1
        lastsymbol = currline(jj-1);
        dofc = ~abs(dofchange(jj+1)-dofchange(jj));
        aofc = aofchange(jj);
        currsblidx = symbolTable(lastsymbol, aofc, dofc);
        if isempty(find(currsblidx))
            sbtable(ii,:)
            aofc
            dofc
            lastsymbol
        end
        curridx = [curridx, find(currsblidx)];
        currstr = [currstr, strlist(currsblidx)];
    end
    minidxarray = curridx;
    for jj = 1:ny-1
        tmpidxarray = circshift(curridx, [0, jj]);
        for kk = 1:ny
            if tmpidxarray(kk) == minidxarray(kk)
                continue;
            end
            if tmpidxarray(kk) < minidxarray(kk)
                minidxarray = tmpidxarray;
                break;
            else
                break; 
            end
        end
    end
    newtable(ii, :) = strlist(minidxarray);
end

    function currsblidx = symbolTable(lastsymbol, aofc, dofc)
        if mod(lastsymbol, 2)
            % last long
            aidx = 1;
        else
            % last short
            aidx = 2;
        end
        if dofc
            % next same
            didx = 3;
        else
            % next opposite
            didx = 4;
        end
        if aofc == 6
            didx = 3;
        end
        currsblidx = ismember(fdtable(:, [aidx didx]), [aofc, 1], 'rows');
    end

end

