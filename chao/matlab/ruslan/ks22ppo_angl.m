% Compute minimum angle of Floquet eigenvectors along a preperiodic orbit.
% ES, 2011-10-09

clear;

% Load periodic orbits
load ks22f90h25;

ph=0; % periodic orbits, only

h=0.25; N=16;

np=5;

global PFLG;  PFLG = -1;

for ipo=1:size(ppo,2),
    
    disp(['iteration ' num2str(ipo)]);
    
    a0=ppo(ipo).a;
    
    [tt, aa, da] = ksfmetd_i(a0, L, h, ppo(ipo).T, np);
    [tt2, aa2, da2] = ksfmetd_i(aa(:,end), L, h, ppo(ipo).T, np);
    % [tt, aa] = ksfmedt(L,tend, a0, h,  np);

%      [x, uu] = ksfm2real(aa, L);
%  
%      fig2 = figure('PaperPosition',[8 8 4 10],'PaperOrient','port',...
%                        'PaperSize',[20.98 29.68],'Position',[400 300 180 450]);
%      pcolor(x,tt,uu'); caxis([-3 3]); shading flat; hold off;

     [ev0, eils] = eig(da(:,end-(2*N-2)+1:end)*da2(:,end-(2*N-2)+1:end));
     eils=diag(eils);
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
    
%     eils=diag(eil);
    
    nflqt=5; % How deep to go into the spectrum;
    unittol=1e-5;

    angl=zeros(size(aa,2),1);

    for i=1:size(aa,2),
        kj=0;
        ev=Jev(da, ev0,i); % use  ev0(:,1:nflqt) ?
        for k=1:nflqt,
            if norm(eils(k)) > 1 && norm(eils(k)-1)>1e-3 % expanding exponent
                if imag(eils(k)) == 0, % real eigenvalue
                    ev(:,k)=ev(:,k)/norm(ev(:,k));
                    for j=k+1:nflqt,
                        ev(:,j)=ev(:,j)/norm(ev(:,j));
                        if norm(eils(j)-1)>1e-3 && imag(eils(j)) == 0 %
                            kj=kj+1;
                            angl(i,kj)=acos(dot(ev(:,j),ev(:,k)));
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
                            angl1=acos(dot(ev(:,k),evComp));
                            if angl1>pi/2., angl1=pi-angl1; end;
                            angl2=acos(dot(ev(:,k),-evComp));
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
                            angl1=acos(dot(ev(:,j),evComp));
                            angl2=acos(dot(ev(:,j),-evComp));
                            angl(i,kj)=min(angl1,angl2);
                        elseif abs(abs(eils(j))-1)>1e-3 && imag(eils(j)) ~= 0
                            exprt= [ipo, ppo(ipo).T, j, eils(j)];
                            save('ks22ppo_angl_probl.dat', 'exprt', '-ascii','-double','-tabs', '-append');
                        end
                    end    
                end
            else
                break; % skip rest of loop when we find first marginal exponent
            end
        end
    end
    
   ppo(ipo).angl =  min(angl);
   exprt= [ppo(ipo).T, min(angl)];
     
    save('ks22ppo_min_angl.dat', 'exprt', '-ascii','-double','-tabs', '-append');
   
   if (mod(ipo,1000)==0), save ks22f90h25angl ppo -append; end; % save in mat file every 1000 orbits.
   
end

save ks22f90h25angl ppo -append;


