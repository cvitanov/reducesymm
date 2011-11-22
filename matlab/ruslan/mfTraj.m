function [aa th] = mfTraj(a,slc)

aa=a;

if nargout < 2,
    for i=1:size(aa,2),
        aa(:,i)=mf(aa(:,i),slc);
    end
else
    th=zeros(size(aa,2),1);
    for i=1:size(aa,2),
        [aa(:,i), th(i)] =mf(aa(:,i),slc);
    end
end



end