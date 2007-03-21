Nd=128;
uexp=[];
for i=1:size(a,2)
    [x,u]=ksfm2real(a(:,i),L,Nd);
    uexp=[uexp -u(1:size(u)-1)/2];
%     iname=['../../production/KS22.0/rpo/' num2str(i) '.dat'];
%    save(iname, 'uexp', '-ASCII');
end 
uexp=transpose(uexp);
save('../../../rpo_ks/davidchack/orbits/rpo.dat','uexp', '-ASCII');