function aa = mfTraj(a,slc)

aa=a;

for i=1:size(aa,2),
    aa(:,i)=mf(aa(:,i),slc);
end

end