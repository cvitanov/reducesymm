function aa = mfinvTraj(a)

aa=zeros(4*size(a,1),size(a,2));
% disp(size(aa));
% disp(size(mfinvSO2(a(:,1))));
for i=1:size(a,2),
    aa(:,i)=mfinvSO2(a(:,i));
end

end