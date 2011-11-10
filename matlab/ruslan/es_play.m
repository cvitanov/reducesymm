clear;
ipo=6; N=16; h=0.25;
load ks22f90h25t100;

format long;

[f,df] = ksfmms3([ppo(ipo).a;ppo(ipo).T;0],L,h); 
disp(['|f|=' num2str(norm(f(1:end-2)))]);
dftmp=df(1:end-2,1:end-2)+eye(2*N-2);
dftmp=dftmp*dftmp;
[vdf,edf] = eig(dftmp); edf = diag(edf);
[sedf, ie] = sort(abs(edf),1,'descend');
ppo(ipo).eig = edf(ie);  ppo(ipo).evec = vdf(:,ie);