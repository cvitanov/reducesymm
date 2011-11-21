function aa = mf1Traj(a)

aa=a;

for i=1:size(aa,2),
    aa(:,i)=mf1(aa(:,i));
end

end