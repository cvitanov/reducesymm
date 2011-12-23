function aa = ksfmpoinc(a, ao, v, dir)

atmp=a-ao*ones(1,size(a,2));
atmp=v'*atmp;
asgn=sign(atmp);
asgn=diff(asgn);
pzero=find(dir*asgn > 1); % we have asign == 1 for point which defines poincare section, asign == 2 when crossing 0
aa=[];
for i=(1:size(pzero,2)),
    aa=[aa, a(:,pzero(i))+(atmp(pzero(i))/(atmp(pzero(i))-atmp(pzero(i)+1)))*(a(:,pzero(i)+1)-a(:,pzero(i)))]; % From linear interpolation
end