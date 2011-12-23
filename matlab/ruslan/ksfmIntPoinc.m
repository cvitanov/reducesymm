function [aPoinc, aout]= ksfmIntPoinc(a0, ao, v, dir, L, h, nhit, tmax)
% ao and v should be in reduced space; a in original space

aPoinc=zeros(size(ao,1),nhit);
nmax=floor(tmax/(4*h));
i=1;
aa=[];
[~, atmp] = ksfmetd2dt(a0,L,h,3,1); % Get enough points for 3rd order interpolation
aa=[aa, atmp];
adum = ksfmPoincInt(aa, ao, v, dir, L, h); 
for j=1:nmax,
    if (size(adum) ~= 0)
        aPoinc(:,i)=mfinvTraj(adum);
        adum=[];
        if i<nhit
            i=i+1;
        else
            if nargout>1
                aout=mfinvTraj(aa);
%                 aout(:,end)=aPoinc(:,i); % "Integrate up to nhit
%                 intersections": does not work properly
            end
            return;
        end
    end
    [~, atmp] = ksfmetd2dt(aa(:,end),L,h,3,1); % Get enough points for 3rd order interpolation
    aa=[aa, atmp(:,2:end)]; % First point of atmp is last point of aa
    adum = ksfmPoincInt(aa(:,end-3:end), ao, v, dir, L, h, aa(:,end-4)); %
end
if (size(adum) ~= 0)
        aPoinc(:,i)=mfinvTraj(adum);
        adum=[];
        if i<nhit
            i=i+1;
        else
            if nargout>1
                aout=mfinvTraj(aa);
%                 aout(:,end)=aPoinc(:,i); % "Integrate up to nhit
%                 intersections": does not work properly
            end
            return;
        end
else
    if nargout>1,
        if (size(adum)==0), aout=[]; end;
    end
end



