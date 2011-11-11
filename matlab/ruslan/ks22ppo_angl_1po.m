% Compute minimum angle of Floquet eigenvectors along a single preperiodic orbit.
% ES, 2011-11-11
    
clear;

% Load periodic orbits
load ks22f90h25;

ph=0; % periodic orbits, only

h=0.25; N=16;

np=1;

global PFLG;  PFLG = -1;

ipo=48;

a0=ppo(ipo).a;

ppo(ipo).angl_prob = false;

[tt, aa, da] = ksfmetd_i(a0, L, h, ppo(ipo).T, np);
a0r=a0;
a0r(1:2:end-1)=-a0r(1:2:end-1);
%     for i=1:2:size(a0r),
%         a0r(i)=-a0r(i);
%     end
[tt2, aa2, da2] = ksfmetd_i(a0r, L, h, ppo(ipo).T, np);%ksfmetd_i(aa(:,end), L, h, ppo(ipo).T, np);
% [tt, aa] = ksfmedt(L,tend, a0, h,  np);

%      [x, uu] = ksfm2real(aa, L);
%  
%      fig2 = figure('PaperPosition',[8 8 4 10],'PaperOrient','port',...
%                        'PaperSize',[20.98 29.68],'Position',[400 300 180 450]);
%      pcolor(x,tt,uu'); caxis([-3 3]); shading flat; hold off;

[ev0, eils] = eig(da2(:,end-(2*N-2)+1:end)*da(:,end-(2*N-2)+1:end));
eils=diag(eils);

% Check accuracy of eigenvalue computation against Ruslan's
% listed values.

if norm(abs(eils)-abs(ppo(ipo).e.^2))>0.1, ppo(ipo).angl_prob = true; end;

%      da2=da(:,end-size(da,1)+1:end);
%      for i=1:2:size(da2,1)
%          for j=1:2:size(da2,2)
%              da2(i,j)=-da2(i,j);
%          end
%      end
%      
%     [ev0, eil] = eig(da(:,end-size(da,1)+1:end)*da2); % Find Floquet vectors/multipliers

err=norm(aa(1:2:end-1,1)+aa(1:2:end-1,end)+1i*(aa(2:2:end,1)-aa(2:2:end,end)));
%     err=norm(aa(1:2:end-1,1)-aa(1:2:end-1,end)+1i*(aa(2:2:end,1)-aa(2:2:end,end)));
disp(['err= ' num2str(err)]);
err2=norm(aa2(:,end)-a0);
disp(['err2= ' num2str(err2)]);
%     eils=diag(eil);

nflqt=5; % How deep to go into the spectrum;
unittol=1e-5;

angl=zeros(size(aa,2),1);
dotp=zeros(size(aa,2),1);

for i=1:size(aa,2),
    kj=0;
    ev=Jev(da, ev0,i); % use  ev0(:,1:nflqt) ?
    for k=1:nflqt,
        if norm(eils(k)) > 1 && norm(eils(k)-1)>1e-3 % expanding exponent
            if imag(eils(k)) == 0, % real eigenvalue
                if imag(ev(1,k)) ~= 0, ppo(ipo).angl_prob = true; end; % mark complex eigenvectors associated to real eigenvalues as problematic
                ev(:,k)=ev(:,k)/norm(ev(:,k));
                for j=k+1:nflqt,
                    ev(:,j)=ev(:,j)/norm(ev(:,j));
                    if norm(abs(eils(j))-1)>1e-3 && imag(eils(j)) == 0 %
                        if imag(ev(1,j)) ~= 0, ppo(ipo).angl_prob = true; end; % mark complex eigenvectors associated to real eigenvalues as problematic
                        kj=kj+1;
                        dp=dot(ev(:,j),ev(:,k));
                        dotp(i,kj)=dp;
                        if abs(dp-1)<2.5*eps, dp =1.; end;
                        if abs(dp+1)<2.5*eps, dp =-1.; end;
                        angl(i,kj)=acos(dp);
                        if angl(i,kj)>pi/2., angl(i,kj)=pi-angl(i,kj); end;
                    elseif abs(abs(eils(j))-1)>1e-3 && imag(eils(j)) ~= 0
                        if eils(j-1)== conj(eils(j)) , % skip step if we have checked cc.
                            continue;
                        end
                        kj=kj+1;
                        ev1=real(ev(:,j))/norm(real(ev(:,j)));
                        ev2=imag(ev(:,j))/norm(imag(ev(:,j)));
                        c =  ( dot(ev2,ev(:,k))-dot(ev1,ev2)*dot(ev1,ev(:,k)) )/...
                        ( dot(ev1,ev(:,k))-dot(ev1,ev2)*dot(ev2,ev(:,k)) );
                        a=1/sqrt(1+c^2+2*c*dot(ev1,ev2));
                        b=c*a;
                        evComp=a*ev1+b*ev2;
                        dp=dot(ev(:,k),evComp);
                        dotp(i,kj)=dp;
                        if abs(dp-1)<2*eps, dp =1.; end;
                        if abs(dp+1)<2*eps, dp =-1.; end;
                        angl(i,kj)=acos(dp);
                        angl1=acos(dp);
                        if angl1>pi/2., angl1=pi-angl1; end;
                        angl2=acos(-dp);
                        if angl2>pi/2., angl2=pi-angl2; end;
                        angl(i,kj)=min(angl1,angl2);
                    end
                end
            else
                if k>1 && eils(k-1) == conj( eils(k) ), % skip step if we have checked cc.
                    continue;
                end
                ev1=real(ev(:,k))/norm(real(ev(:,k)));
                ev2=imag(ev(:,k))/norm(imag(ev(:,k)));
                for j=k+2:nflqt,
                    c = ( dot(ev2,ev(:,j))-dot(ev1,ev2)*dot(ev1,ev(:,j)) )/...
                        ( dot(ev1,ev(:,j))-dot(ev1,ev2)*dot(ev2,ev(:,j)) );
                    a=1/sqrt(1+c^2+2*c*dot(ev1,ev2));
                    b=c*a;
                    evComp=a*ev1+b*ev2;
                    if abs(abs(eils(j))-1)>1e-3 && imag(eils(j)) == 0, % then we would have to check angle of subspaces
                        ev(:,j)=ev(:,j)/norm(ev(:,j));
                        kj=kj+1;
                        dp=dot(ev(:,j),evComp);
                        dotp(i,kj)=dp;                                     
                        if abs(dp-1)<2*eps, dp =1.; end;
                        if abs(dp+1)<2*eps, dp =-1.; end;
                        angl(i,kj)=acos(dp);
                        angl1=acos(dp);
                        if angl1>pi/2., angl1=pi-angl1; end;
                        angl2=acos(-dp);
                        if angl2>pi/2., angl2=pi-angl2; end;
                        angl(i,kj)=min(angl1,angl2);
                        if imag(angl(i,kj)) ~= 0, ppo(ipo).angl_prob = true; end;          
                    elseif abs(abs(eils(j))-1)>1e-3 && imag(eils(j)) ~= 0
%                              exprt= [ipo, ppo(ipo).T, j, eils(j)];
                        ppo(ipo).angl_prob = true;
%                             save('ks22ppo_angl_probl.dat', 'exprt', '-ascii','-double','-tabs', '-append');
                    end
                end    
            end
        else
            break; % skip rest of loop when we find first marginal exponent
        end
        if (kj>0 && imag(angl(i,kj)) ~= 0), ppo(ipo).angl_prob = true; end;   % mark complex valued angles as problematic
    end
end

     
save(['ks22ppo_angl_', num2str(48), '.dat'], 'angl', '-ascii','-double','-tabs');
