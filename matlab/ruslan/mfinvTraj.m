function aa = mfinvTraj(a)

aa=a;

for i=1:size(aa,2),
    aa(:,i)=mfinvSO2(aa(:,i));
end

end