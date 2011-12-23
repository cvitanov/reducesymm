function rm = manifRep2rm(a, ao, aman, s, io, cut)
% Return map for "first fold" of unstable manifold.
% a: points on the manifold (reduced space) with unknown s.
% ao: point of reference for s [s(ao)=0]
% aman: part of the manifold for which s has been computed
% s: tabulated values of s
% io: index of first point of the previous iteration of the manifold
% cut: distance cutoff

rm =[];% zeros(size(a,2),2);
for i=1:size(a,2),
    stmp=p2s(a(:,i), ao, aman, s, cut);
    if not(isnan(stmp)),
        disp(i);
        rm = [rm, [s(io+1+i) stmp]']; % 
    end    
end

rm=rm';

end

